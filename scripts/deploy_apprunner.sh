#!/usr/bin/env bash
# Put the question server on AWS App Runner, from the committed code, or update it to a new commit.
#
# Builds the x86-64 image from `git archive HEAD` and the release archive (the image checks every file
# of the bundle against the committed digests), pushes it to an ECR repository, and creates or updates
# one App Runner service: 1 vCPU and 2 GB, never more than one instance, health-checked at /health.
# App Runner gives it an HTTPS address; the MCP endpoint is that address with /mcp.
#
#   scripts/deploy_apprunner.sh            # needs the AWS CLI signed in, docker with buildx, and
#                                          # data/release/pharos-gap-kg-bundle-pharos319.tar.gz
# Cost while running (us-east-1, AWS price list, 27 Oct 2026): 2 GB kept on at $0.007 per GB-hour,
# plus $0.064 per vCPU-hour while answering. `aws apprunner pause-service` stops the charge for compute.
set -euo pipefail
NAME=pharos-gap-kg
REGION="${AWS_REGION:-us-east-1}"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
ARCHIVE=data/release/pharos-gap-kg-bundle-pharos319.tar.gz
if [ -n "$(git status --porcelain -- src pyproject.toml constraints.txt Dockerfile .dockerignore docs/data/bundle_archive_pharos319.json)" ]; then
  echo "REFUSED: what the image is built from has uncommitted changes; the image is the commit" >&2; exit 1
fi
[ -f "$ARCHIVE" ] || { echo "missing $ARCHIVE; fetch it:  python -m pharos_gap_kg.release get pharos319" >&2; exit 1; }
TAG="$(git rev-parse --short HEAD)"
ACCOUNT="$(aws sts get-caller-identity --query Account --output text)"
REGISTRY="$ACCOUNT.dkr.ecr.$REGION.amazonaws.com"
IMAGE="$REGISTRY/$NAME:$TAG"

echo "1. the image, from $TAG"
SRC="$(mktemp -d)"; trap 'rm -rf "$SRC"' EXIT
git archive HEAD | tar -x -C "$SRC"
mkdir -p "$SRC/data/release" && cp "$ARCHIVE" "$SRC/data/release/"
docker buildx build --platform linux/amd64 --provenance=false --sbom=false -t "$IMAGE" --load "$SRC"

echo "2. the registry"
aws ecr describe-repositories --region "$REGION" --repository-names "$NAME" >/dev/null 2>&1 \
  || aws ecr create-repository --region "$REGION" --repository-name "$NAME" --image-scanning-configuration scanOnPush=true >/dev/null
aws ecr get-login-password --region "$REGION" | docker login --username AWS --password-stdin "$REGISTRY" >/dev/null
docker push -q "$IMAGE"
docker logout "$REGISTRY" >/dev/null   # the login token is written to ~/.docker/config.json unencrypted

echo "3. the role App Runner reads the registry with"
ROLE="AppRunnerECRAccessRole-$NAME"
if ! aws iam get-role --role-name "$ROLE" >/dev/null 2>&1; then
  aws iam create-role --role-name "$ROLE" --assume-role-policy-document \
    '{"Version":"2012-10-17","Statement":[{"Effect":"Allow","Principal":{"Service":"build.apprunner.amazonaws.com"},"Action":"sts:AssumeRole"}]}' >/dev/null
  aws iam attach-role-policy --role-name "$ROLE" --policy-arn arn:aws:iam::aws:policy/service-role/AWSAppRunnerServicePolicyForECRAccess
  sleep 15   # a new role takes a few seconds to be usable
fi
ROLE_ARN="$(aws iam get-role --role-name "$ROLE" --query Role.Arn --output text)"

echo "4. at most one instance"
SCALE="$(aws apprunner list-auto-scaling-configurations --region "$REGION" --auto-scaling-configuration-name "$NAME-one" \
  --query 'AutoScalingConfigurationSummaryList[0].AutoScalingConfigurationArn' --output text)"
if [ "$SCALE" = "None" ] || [ -z "$SCALE" ]; then
  SCALE="$(aws apprunner create-auto-scaling-configuration --region "$REGION" --auto-scaling-configuration-name "$NAME-one" \
    --min-size 1 --max-size 1 --max-concurrency 50 --query AutoScalingConfiguration.AutoScalingConfigurationArn --output text)"
fi

echo "5. the service"
SOURCE="{\"ImageRepository\":{\"ImageIdentifier\":\"$IMAGE\",\"ImageRepositoryType\":\"ECR\",\"ImageConfiguration\":{\"Port\":\"8080\"}},
  \"AutoDeploymentsEnabled\":false,\"AuthenticationConfiguration\":{\"AccessRoleArn\":\"$ROLE_ARN\"}}"
SERVICE="$(aws apprunner list-services --region "$REGION" --query "ServiceSummaryList[?ServiceName=='$NAME'].ServiceArn | [0]" --output text)"
if [ "$SERVICE" = "None" ] || [ -z "$SERVICE" ]; then
  SERVICE="$(aws apprunner create-service --region "$REGION" --service-name "$NAME" --source-configuration "$SOURCE" \
    --instance-configuration Cpu=1024,Memory=2048 \
    --health-check-configuration Protocol=HTTP,Path=/health,Interval=10,Timeout=5,HealthyThreshold=1,UnhealthyThreshold=5 \
    --auto-scaling-configuration-arn "$SCALE" --query Service.ServiceArn --output text)"
else
  aws apprunner update-service --region "$REGION" --service-arn "$SERVICE" --source-configuration "$SOURCE" >/dev/null
fi

echo "6. waiting for it to run"
for _ in $(seq 1 90); do
  STATUS="$(aws apprunner describe-service --region "$REGION" --service-arn "$SERVICE" --query Service.Status --output text)"
  case "$STATUS" in
    RUNNING) break ;;
    CREATE_FAILED|DELETED|PAUSED) echo "the service is $STATUS; see its event log in the App Runner console" >&2; exit 1 ;;
  esac
  sleep 20
done
URL="$(aws apprunner describe-service --region "$REGION" --service-arn "$SERVICE" --query Service.ServiceUrl --output text)"
echo "running $TAG at https://$URL/mcp  (health: https://$URL/health)"
