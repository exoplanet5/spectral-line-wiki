"""Print id/title/journal-ref/doi/abstract from an arXiv API Atom XML file (untrusted input; run with python -I)."""
import sys
import xml.etree.ElementTree as ET

ns = {'a': 'http://www.w3.org/2005/Atom', 'x': 'http://arxiv.org/schemas/atom'}
for path in sys.argv[1:]:
    root = ET.parse(path).getroot()
    for e in root.findall('a:entry', ns):
        def g(tag):
            el = e.find(tag, ns)
            return ' '.join(el.text.split()) if el is not None and el.text else None
        print('ID   ', g('a:id'))
        print('TITLE', g('a:title'))
        print('AUTH ', ', '.join(' '.join(a.find('a:name', ns).text.split()) for a in e.findall('a:author', ns)[:6]))
        print('JREF ', g('x:journal_ref'))
        print('DOI  ', g('x:doi'))
        print('ABS  ', g('a:summary'))
        print()
