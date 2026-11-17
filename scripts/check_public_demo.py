#!/usr/bin/env python3
"""Verify the submission's public pages, functional APIs, MCP and downloads."""
import argparse
import asyncio
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import time
import urllib.error
import urllib.request
from urllib.parse import quote, urljoin, urlparse


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--base', default='https://pm4ie6tprh.us-east-1.awsapprunner.com')
    p.add_argument('--commit', required=True)
    p.add_argument('--tag', default='submission-2026-11-17-r3')
    p.add_argument('--output', required=True)
    a = p.parse_args()
    results = []
    def get(url):
        req = urllib.request.Request(url)
        for attempt in range(3):
            try:
                with urllib.request.urlopen(req, timeout=60) as r:
                    data = r.read()
                    assert r.status == 200, (url, r.status)
                    results.append({'url': url, 'status': r.status, 'bytes': len(data)})
                    return data
            except urllib.error.HTTPError as e:
                results.append({'url': url, 'status': e.code, 'attempt': attempt + 1})
                if e.code not in (429, 502, 503, 504) or attempt == 2:
                    raise
                time.sleep(2 * (attempt + 1))
    def check():
        h = json.loads(get(a.base + '/health'))
        assert h['status'] == 'ok' and h['commit'] == a.commit, h
        for route in ['/', '/add?run=D']:
            html = get(a.base + route).decode()
            assert '<html' in html.lower()
            for ref in sorted(set(re.findall(r'(?:src|href)=["\']([^"\']+)["\']', html))):
                if ref.endswith(('.css', '.js')):
                    u = urljoin(a.base, ref)
                    assert urlparse(u).netloc == urlparse(a.base).netloc
                    get(u)
        walk = json.loads(get(a.base + '/walkthrough.json'))
        assert walk['studies'], 'walkthrough has no study examples'
        if a.tag in {'submission-2026-11-17-r2', 'submission-2026-11-17-r3'}:
            study = next(s for s in walk['studies'] if s['id'] == 'D')
            assert study['join']['targets'] == 500 and '499 genes' in study['reference_note'], study['id']
            why = json.loads(get(a.base + '/api/why?target=ACHE&trait=Alzheimer%20disease'))
            assert 'donepezil' in why['answer'].casefold(), 'missing drug path in the featured example'
        found = json.loads(get(a.base + '/api/find?text=Parkinson%20disease'))
        assert found.get('matches'), found
        gaps = json.loads(get(a.base + '/api/gaps?disease=Parkinson%20disease'))
        assert gaps.get('found') is True, gaps
        if a.tag == 'submission-2026-11-17-r3':
            assert (gaps['holds']['tied_genes'], gaps['counts']['no_genetic_support']) == (345, 97)
            assert '55,000' in gaps['answer']
            diabetes = json.loads(get(a.base + '/api/gaps?disease=type%202%20diabetes%20mellitus'))
            assert diabetes['holds']['tied_genes'] == 2549 and diabetes['counts']['no_genetic_support'] == 115
            supporting = json.loads(get(a.base + '/api/subgraph?question=gaps&disease=MONDO%3A0005180'))
            assert supporting['nodes'] and supporting['edges']
        status = json.loads(get(a.base + '/api/status'))
        assert status, 'empty graph status'
        async def mcp():
            from mcp import ClientSession
            from mcp.client.streamable_http import streamable_http_client
            async with streamable_http_client(a.base + '/mcp') as streams:
                async with ClientSession(streams[0], streams[1]) as s:
                    await s.initialize()
                    tools = [t.name for t in (await s.list_tools()).tools]
                    assert 'gaps' in tools and not any(t.startswith('add_') for t in tools), tools
                    r = await s.call_tool('gaps', {'disease': 'Parkinson disease'})
                    assert not getattr(r, "is_error", getattr(r, "isError", True)), r
                    body = getattr(r, 'structured_content', None) or getattr(r, 'structuredContent', None)
                    if body is None:
                        body = json.loads(r.content[0].text)
                    assert body.get('found') is True, body
                    results.append({'url': a.base + '/mcp', 'initialize': 'passed', 'tools': tools, 'gaps_call': 'passed'})
        asyncio.run(mcp())
        api = 'https://api.github.com/repos/iandriver/pharos-gap-kg'
        tag = quote(a.tag, safe='')
        repo = 'https://github.com/iandriver/pharos-gap-kg'
        commit = json.loads(get(api + '/commits/' + tag))
        assert commit['sha'] == a.commit, commit['sha']
        get(repo)
        get(f'{repo}/tree/{tag}')
        for name in ['SUBMISSION_20261117.md', 'PROGRESS_20261027.md', 'PROTOCOL.md', 'IMPORT_WORKFLOW.md', 'VIDEO.md']:
            get(f'{repo}/blob/{tag}/docs/{name}')
        get(f'{repo}/blob/{tag}/findings/REPRODUCE.md')
        if a.tag == 'submission-2026-11-17-r3':
            assets = json.loads(get(api + '/releases/tags/' + tag))['assets']
            records = {x['name']: x for x in assets}
            for name in ['pharos-gap-kg-final-submission.docx', 'pharos-gap-kg-gaps-demo.mp4']:
                asset = records[name]
                content = get(asset['browser_download_url'])
                digest = hashlib.sha256(content).hexdigest()
                assert asset['digest'] == 'sha256:' + digest, name
                results[-1]['sha256'] = digest
        release = json.loads(get(api + '/releases/tags/pharos319-bundle'))
        assets = {x['name']: x for x in release['assets']}
        checks = {'pharos-gap-kg-graph-pharos319.tar.gz': '4891f7c80671a5771348817af16a8c9a8bc2652e4792db14d5221e7474d8c72c',
                  'pharos-gap-kg-bundle-pharos319.tar.gz': 'a6b1554c9cac320836652f1f198920c3013491ccc8081fbfd38ba112d8b0a3d3'}
        for name, wanted in checks.items():
            blob = get(assets[name]['browser_download_url'])
            observed = hashlib.sha256(blob).hexdigest()
            assert observed == wanted, (name, observed, wanted)
            results[-1]['sha256'] = observed
        checksum = get(assets['pharos-gap-kg-graph-pharos319.tar.gz.sha256']['browser_download_url']).decode()
        assert checks['pharos-gap-kg-graph-pharos319.tar.gz'] in checksum
    passed, error = False, None
    try:
        check()
        passed = True
    except Exception as e:
        error = f'{type(e).__name__}: {e}'
    out = {'checked_at': datetime.now(timezone.utc).isoformat(), 'commit': a.commit, 'tag': a.tag, 'passed': passed,
           'checks': results, 'error': error}
    Path(a.output).write_text(json.dumps(out, indent=2) + '\n')
    print(json.dumps({'passed': passed, 'checks': len(results), 'error': error, 'output': a.output}))
    return 0 if passed else 1


if __name__ == '__main__':
    raise SystemExit(main())
