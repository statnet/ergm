import re
from operator import add
from collections import defaultdict
import sys

HEAD_RE = re.compile("^(?P<t>[@+ -])(@ *-(?P<s>[0-9]+))?")

class LineMapper:
    def __init__(self, patch):
        matches = [HEAD_RE.match(line) for line in patch.split('\n')]
        starts = [int(m['s']) for m in matches if m['t'] == '@']
        lines = ''.join([m['t'] for m in matches[1:]]).split('@')
        lens = [len(l) for l in lines]
        ends = list(map(add, starts, lens))

        self.lmap = [0] * max(ends)

        pos1 = 0
        pos2 = 0

        for ch in range(len(starts)):
            gap = starts[ch] - pos1 - 1
            self.lmap[pos1:(pos1 + gap)] = range(pos2 + 1, pos2 + gap + 1)

            pos1 = pos1 + gap
            pos2 = pos2 + gap

            for d in lines[ch]:
                if d in [' ', '-']: pos1 += 1
                if d in [' ', '+']: pos2 += 1
                if d == ' ': self.lmap[pos1 - 1] = pos2

        gap = len(self.lmap) - pos1
        self.lmap[pos1:(pos1 + gap)] = range(pos2 + 1, pos2 + gap + 1)

        self.shift = pos2 - pos1

    def __getitem__(self, key):
        if key == 0: return 0
        elif key <= len(self.lmap): return self.lmap[key - 1]
        else: return key + self.shift

def chunk(l, each):
    for i in range(0, len(l), each):
        yield l[i:i+each]

def split_lines(lines, prefix, each = 10): # 10 = GitHub max per step.
    for i, l in enumerate(chunk(lines, each)):
        open(prefix + str(i + 1) + '.txt', 'w').writelines(l)

def split_by_file(lines, file_regex):
    o = defaultdict(list)
    for l in lines:
        f = file_regex.match(l)['file']
        o[f].append(l)

    return dict(o)

def split_patch_by_file(patch):
    patch = patch.split('\n')
    o = {}
    f = None

    for l in patch:
        if l.startswith('diff'):
            parts = l.strip().split()
            # The file name is usually after 'b/' (the new file)
            for part in parts:
                if part.startswith('b/'):
                    f = part[2:]
        else:
            o[f] += l + '\n'

    return dict(o)

def uniq_messages(old, new):
    return [msg for msg in new if msg not in set(old)]

def remap_file_messages(old, patch, line_regex, line_replace):
    m = LineMapper(patch)
    ol = [int(line_regex.match(l)['line']) for l in old]
    nl = [m[l] for l in ol]

    return list(map(line_replace, old, ol, nl))

def filter_messages(old, new, patch_dict, message_regex, line_replace):
    old = split_by_file(old, message_regex)
    new = split_by_file(new, message_regex)

    o = []
    for f in new.keys() & patch_dict.keys():
        oldnew = remap_file_messages(old.get(f, []), patch_dict[f], message_regex, line_replace)
        o.extend(uniq_messages(oldnew, new[f]))

    return o

