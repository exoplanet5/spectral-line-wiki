#!/usr/bin/env python
"""Headless-Chrome smoke test of index.html.
usage: python tests/headless.py [screenshot.png [width height [#hash]]]
Prints the JSON self-test report written by the page (?selftest) and optionally a screenshot."""
import html, json, os, re, subprocess, sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
CH = '/Applications/Google Chrome.app/Contents/MacOS/Google Chrome'
args = sys.argv[1:]
PAGE = os.environ.get('SLW_PAGE', os.path.join(ROOT, 'index.html'))
url = 'file://' + PAGE + '?selftest' + os.environ.get('SLW_QS', '') + ('&scroll=' + args[3] if len(args) > 3 else '')
dom = subprocess.run([CH, '--headless=new', '--disable-gpu', '--virtual-time-budget=9000', '--dump-dom', url],
                     capture_output=True, text=True, timeout=180).stdout
m = re.search(r'<pre id="selftest">(.*?)</pre>', dom, re.S)
print(json.dumps(json.loads(html.unescape(m.group(1))), indent=1) if m else 'NO SELFTEST OUTPUT')
if args:
    w, h = (args[1], args[2]) if len(args) > 2 else ('1400', '2400')
    subprocess.run([CH, '--headless=new', '--disable-gpu', '--hide-scrollbars', '--virtual-time-budget=9000',
                    f'--window-size={w},{h}', f'--screenshot={args[0]}', url], capture_output=True, timeout=180)
    print('screenshot', args[0])
