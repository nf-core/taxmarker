#!/usr/bin/env python3
"""Make the dark metro map from the light one: python3 make_dark_metro_map.py light.svg dark.svg

The dark map has the page of the light one, a dark background and inverted ink.
Colours that are not in the mapping below are kept; an element without a fill is light.
"""
import re, sys
src, dst = sys.argv[1:3]
s = open(src).read()
BG, FG, LABEL = '#2b2b2b', '#f2f2f2', '#e6e6e6'
def prop(name, mapping):
    global s
    for old, new in mapping.items():
        s, n = re.subn(rf'(?i)((?<![-\w]){name}\s*[:=]\s*"?)#{old[1:]}\b', lambda m: m.group(1) + new, s)
        print(f'{name} {old} -> {new}: {n}')
prop('fill', {'#ffffff': BG, '#000000': FG, '#333333': LABEL, '#0b0908': '#fafafa', '#3f2b29': '#764f4b'})
prop('stroke', {'#000000': FG})
s, n = re.subn(r'(<g\s+inkscape:label="Layer 1"\s+inkscape:groupmode="layer")', r'\1\n     fill="' + FG + '"', s, count=1)
print('layer default fill:', n)
x, y, w, h = re.search(r'<svg\b[^>]*?viewBox="([-\d.]+) ([-\d.]+) ([\d.]+) ([\d.]+)"', s, re.S).groups()
s = s.replace('pagecolor="#ffffff"', f'pagecolor="{BG}"')
s = re.sub(r'sodipodi:docname="[^"]*"', 'sodipodi:docname="' + dst.split('/')[-1] + '"', s)
rect = f'<rect id="background" x="{x}" y="{y}" width="{w}" height="{h}" fill="{BG}"/>\n  '
i = s.index('<g', s.index('</defs>'))
open(dst, 'w').write(s[:i] + rect + s[i:])
