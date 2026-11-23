"""Minimal HTML -> text converter for inspecting downloaded (untrusted) pages.
Usage: python -I html2text.py <file.html>
"""
import sys, re, html

def convert(t):
    t = re.sub(r'<script.*?</script>', '', t, flags=re.S | re.I)
    t = re.sub(r'<style.*?</style>', '', t, flags=re.S | re.I)
    t = re.sub(r'<img[^>]*alt="([^"]*)"[^>]*>', r' [IMG alt=\1] ', t, flags=re.I)
    t = re.sub(r'<(br|p|div|tr|li|h\d|table)[^>]*>', '\n', t, flags=re.I)
    t = re.sub(r'<(td|th)[^>]*>', ' | ', t, flags=re.I)
    t = re.sub(r'<sup>', '^', t, flags=re.I)
    t = re.sub(r'<sub>', '_', t, flags=re.I)
    t = re.sub(r'<[^>]+>', '', t)
    t = html.unescape(t)
    t = re.sub(r'[ \t]+', ' ', t)
    t = re.sub(r'\n\s*\n+', '\n', t)
    return t

if __name__ == '__main__':
    print(convert(open(sys.argv[1], encoding='utf-8', errors='replace').read()))
