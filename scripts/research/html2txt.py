import re, sys, html
t = open(sys.argv[1], encoding='utf-8', errors='replace').read()
t = re.sub(r'<script.*?</script>', ' ', t, flags=re.S)
t = re.sub(r'<style.*?</style>', ' ', t, flags=re.S)
t = re.sub(r'<[^>]+>', ' ', t); t = html.unescape(t)
t = re.sub(r'[ \t]+', ' ', t)
lines = [l.strip() for l in t.splitlines() if l.strip()]
a = int(sys.argv[2]) if len(sys.argv) > 2 else 0
b = int(sys.argv[3]) if len(sys.argv) > 3 else len(lines)
print('\n'.join(lines[a:b]))
