#!/usr/bin/env python3
"""Divide a constrained tree search into independent parts, and graft the part trees back together.

split: the parts are chosen bottom-up on the partition tree (default: the constraint tree). A node holding more than
--max tips closes its largest children into parts of their own; a closed part is represented in its parent by --reps
tips. When closing children is not enough (many small children of a flat node), consecutive children are grouped.
The last part, the root, holds the representatives of everything below it. Every part is a set of original tips plus
--outgroups tips from neighbouring parts, with the constraint tree pruned to them and the part (and the
representatives of its child parts) forced to be clades. A part whose constraint tree is binary has nothing to
search: its tree is written to <outdir>/resolved/ instead of an alignment and a constraint.

graft: takes the tree searched for each part, cuts out the part (rooted on its outgroup) and puts it in place of its
representatives in the tree of the parent part, starting from the root. The result must keep every clade of the
constraint tree.
"""
import argparse
import json
import os
import random
import sys
import threading
import traceback

__version__ = "1.0.0"


class Node:
    __slots__ = ("name", "length", "kids")

    def __init__(self, name=None, length=None, kids=None):
        self.name, self.length, self.kids = name, length, kids or []


def _split_label(label):
    if label.startswith("'"):
        end = label.rindex("'")
        rest = label[end + 1:]
        return label[1:end].replace("''", "'"), float(rest[1:]) if rest.startswith(":") and len(rest) > 1 else None
    name, _, length = label.partition(":")
    return name or None, float(length) if length else None


def _tokens(text):
    out, cur, quoted, comment = [], "", False, False
    for ch in text:
        if comment:
            comment = ch != "]"
        elif quoted:
            cur += ch
            quoted = ch != "'"
        elif ch == "'":
            cur += ch
            quoted = True
        elif ch == "[":
            comment = True
        elif ch.isspace():
            continue
        elif ch in "(),;":
            if cur:
                out.append(cur)
                cur = ""
            out.append(ch)
        else:
            cur += ch
    return out


def parse(text):
    tok, pos = _tokens(text), 0
    if not tok:
        sys.exit("No tree found in the input")

    def node():
        nonlocal pos
        n = Node()
        if tok[pos] == "(":
            pos += 1
            while True:
                n.kids.append(node())
                t = tok[pos]
                pos += 1
                if t == ")":
                    break
            if pos < len(tok) and tok[pos] not in "(),;":
                n.name, n.length = _split_label(tok[pos])
                pos += 1
        else:
            n.name, n.length = _split_label(tok[pos])
            pos += 1
        return n

    return node()


def _quote(name):
    return "'%s'" % name.replace("'", "''") if any(c in "(),:;[]' \t" for c in name) else name


def write(n, lengths=False):
    s = _quote(n.name) if not n.kids else "(" + ",".join(write(k, lengths) for k in n.kids) + ")"
    if lengths and n.length is not None:
        s += ":%.8f" % n.length
    return s


def post_order(root):
    out, stack = [], [(root, False)]
    while stack:
        x, done = stack.pop()
        if done or not x.kids:
            out.append(x)
        else:
            stack.append((x, True))
            stack.extend((k, False) for k in reversed(x.kids))
    return out


def leaves(root):
    return [x.name for x in post_order(root) if not x.kids]


def prune(n, keep):
    if not n.kids:
        return Node(n.name) if n.name in keep else None
    kids = [k for k in (prune(c, keep) for c in n.kids) if k]
    if not kids:
        return None
    return kids[0] if len(kids) == 1 else Node(kids=kids)


def counts(root, want):
    """Per node id: (number of leaves in `want`, number of leaves)."""
    c = {}
    for x in post_order(root):
        if not x.kids:
            c[id(x)] = (int(x.name in want), 1)
        else:
            c[id(x)] = (sum(c[id(k)][0] for k in x.kids), sum(c[id(k)][1] for k in x.kids))
    return c


def ensure_clade(root, want):
    """Make `want` a clade of `root` by grouping whole children of the lowest node holding it."""
    want = set(want)
    if len(want) < 2:
        return
    c = counts(root, want)
    node = root
    while True:
        nxt = [k for k in node.kids if c[id(k)][0] == len(want)]
        if not nxt:
            break
        node = nxt[0]
    if c[id(node)][1] == len(want):
        return
    inside = [k for k in node.kids if c[id(k)][0]]
    if any(c[id(k)][0] != c[id(k)][1] for k in inside):
        sys.exit("The constraint tree does not allow %s... to form a clade; a partition tree must not conflict with the constraint tree" % sorted(want)[:3])
    node.kids = [k for k in node.kids if not c[id(k)][0]] + [Node(kids=inside)]


def collapse(n):
    n.kids = [collapse(k) for k in n.kids]
    if len(n.kids) == 1:
        child = n.kids[0]
        if n.length is not None and child.length is not None:
            child.length += n.length
        return child
    return n


def reroot_above(root, tip):
    """Copy of the tree rooted on the branch above leaf `tip`."""
    adj = {}
    for x in post_order(root):
        adj.setdefault(id(x), (x, []))
        for k in x.kids:
            adj.setdefault(id(k), (k, []))
            adj[id(x)][1].append((k, k.length))
            adj[id(k)][1].append((x, k.length))
    leaf = next(x for x, _ in adj.values() if not x.kids and x.name == tip)
    (par, length), = adj[id(leaf)][1]
    half = length / 2 if length is not None else None
    top = Node(par.name, half)
    stack = [(par, leaf, top)]
    while stack:
        x, frm, new = stack.pop()
        for y, ln in adj[id(x)][1]:
            if y is not frm:
                child = Node(y.name, ln)
                new.kids.append(child)
                stack.append((y, x, child))
    return Node(kids=[Node(leaf.name, half), collapse(top)])


def find_clade(tree, tips, outside):
    """(tree, node): a tree with a node whose leaves are exactly `tips`, rerooted above `outside` if needed."""
    want = set(tips)
    for t in (tree, None):
        if t is None:
            t = reroot_above(tree, outside)
        c = counts(t, want)
        for x in post_order(t):
            if c[id(x)][0] == c[id(x)][1] == len(want):
                return t, x
    sys.exit("The clade %s... is not monophyletic in the searched tree" % sorted(want)[:3])


def read_fasta(path):
    records, name, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    records[name] = "".join(buf)
                if not line[1:].split():
                    sys.exit("A sequence in the alignment has no name")
                name, buf = line[1:].split()[0], [line]
                if name in records:
                    sys.exit("Duplicate sequence name %s in the alignment" % name)
            elif name is not None:
                buf.append(line if line.endswith("\n") else line + "\n")
    if name is not None:
        records[name] = "".join(buf)
    return records


def is_resolved(root):
    """True when the tree is binary (an unrooted tree has a root of two or three children): RAxML-NG rejects it as a constraint."""
    return len(root.kids) <= 3 and all(len(x.kids) == 2 for x in post_order(root) if x.kids and x is not root)


def split(a):
    if a.min < 2 or a.reps < 1 or a.outgroups < 1:
        sys.exit("--min must be at least 2, --reps and --outgroups at least 1")
    constraint = parse(open(a.constraint).read())
    if len(set(leaves(constraint))) != len(leaves(constraint)):
        sys.exit("The constraint tree has duplicate tip names")
    ptree = parse(open(a.partition_tree).read()) if a.partition_tree else constraint
    if set(leaves(ptree)) != set(leaves(constraint)):
        sys.exit("The partition tree and the constraint tree have different tips")
    alignment = read_fasta(a.alignment)
    missing = set(leaves(constraint)) - set(alignment)
    if missing:
        sys.exit("Tips without a sequence in the alignment: %s" % sorted(missing)[:3])

    parent, residual, jobs = {}, {}, []

    def newjob(tips, at, others, root=False):
        reps = list(tips) if len(tips) <= a.reps else [tips[int(i * len(tips) / a.reps)] for i in range(a.reps)]
        jobs.append({"tips": tips, "reps": reps, "at": at, "others": others, "root": root})
        return reps

    def process(v):
        if not v.kids:
            residual[id(v)] = [v.name]
            return residual[id(v)]
        items = []
        for k in v.kids:
            parent[id(k)] = v
            items.append(process(k))
        total = sum(len(i) for i in items)
        for i in sorted(range(len(items)), key=lambda i: -len(items[i])):
            if total <= a.max:
                break
            if len(items[i]) >= a.min:
                reps = newjob(items[i], v, [items[j] for j in range(len(items)) if j != i])
                total -= len(items[i]) - len(reps)
                items[i] = reps
        while total > a.max:
            groups, g, size = [], [], 0
            for i, it in enumerate(items):
                if g and size + len(it) > a.max:
                    groups.append(g)
                    g, size = [], 0
                g.append(i)
                size += len(it)
            groups = [g for g in groups + [g] if len(g) > 1]
            if not groups:
                # Last resort: an item below --min is closed too, as nothing else can shrink the node.
                i = max(range(len(items)), key=lambda i: len(items[i]))
                if len(items[i]) <= a.reps:
                    sys.exit("A node cannot be reduced below --max %d; raise --max or lower --reps" % a.max)
                reps = newjob(items[i], v, [items[j] for j in range(len(items)) if j != i])
                total -= len(items[i]) - len(reps)
                items[i] = reps
                continue
            merged = {}
            for g in groups:
                tips = [t for i in g for t in items[i]]
                reps = newjob(tips, v, [items[j] for j in range(len(items)) if j not in g])
                total -= len(tips) - len(reps)
                merged[g[0]] = reps
            drop = {i for g in groups for i in g[1:]}
            items = [merged.get(i, it) for i, it in enumerate(items) if i not in drop]
        residual[id(v)] = [t for i in items for t in i]
        return residual[id(v)]

    top = process(ptree)
    newjob(top, ptree, [], root=True)

    def outgroups(j):
        if j["root"]:
            return []
        og = [o[0] for o in j["others"]][:a.outgroups]
        v = j["at"]
        while len(og) < a.outgroups and id(v) in parent:
            p = parent[id(v)]
            og += [residual[id(k)][0] for k in p.kids if k is not v][:a.outgroups - len(og)]
            v = p
        return og

    os.makedirs(a.outdir, exist_ok=True)
    plan = {"jobs": []}
    for i, j in enumerate(jobs):
        og = outgroups(j)
        tips, tipset = j["tips"], set(j["tips"])
        tree = prune(constraint, tipset | set(og))
        for clade in ([tips] if og else []) + [o["reps"] for o in jobs[:i] if set(o["reps"]) <= tipset]:
            ensure_clade(tree, clade)
        job = "job_%03d" % i
        done = is_resolved(tree)
        if done:
            os.makedirs("%s/resolved" % a.outdir, exist_ok=True)
            with open("%s/resolved/%s.resolved.nwk" % (a.outdir, job), "w") as fo:
                fo.write(write(tree) + ";\n")
        else:
            with open("%s/%s.fasta" % (a.outdir, job), "w") as fo:
                fo.writelines(alignment[t] for t in tips + og)
            with open("%s/%s.nwk" % (a.outdir, job), "w") as fo:
                fo.write(write(tree) + ";\n")
        plan["jobs"].append({"id": job, "n": len(tips), "reps": j["reps"], "outgroups": og, "root": j["root"], "resolved": done})
    with open(a.plan, "w") as fo:
        json.dump(plan, fo, indent=1)
    sizes = sorted((j["n"] for j in plan["jobs"]), reverse=True)
    print("%d parts, largest %s, smallest %d, root %d tips" % (len(jobs), sizes[:5], sizes[-1], len(top)))


def check(constraint, tree):
    """Number of clades of the constraint tree that are not monophyletic in `tree`."""
    tips = leaves(constraint)
    if sorted(tips) != sorted(leaves(tree)):
        sys.exit("The grafted tree and the constraint tree have different tips")
    rnd = random.Random(1)
    h = {t: rnd.getrandbits(64) for t in tips}
    mask = (1 << 64) - 1

    def hashes(root):
        memo, size, out = {}, {}, {}
        for x in post_order(root):
            if not x.kids:
                memo[id(x)], size[id(x)] = h[x.name], 1
            else:
                memo[id(x)] = sum(memo[id(k)] for k in x.kids) & mask
                size[id(x)] = sum(size[id(k)] for k in x.kids)
                out[memo[id(x)]] = size[id(x)]
        return out

    have, total = hashes(tree), sum(h.values()) & mask
    # A clade of one tip or of all tips but one says nothing in an unrooted tree.
    return len([x for x, n in hashes(constraint).items() if 2 <= n <= len(tips) - 2 and x not in have and (total - x) & mask not in have])


def graft(a):
    plan = json.load(open(a.plan))
    result = {}
    for j in plan["jobs"]:
        found = [r for r in a.results if os.path.basename(r).startswith(j["id"] + ".")]
        if len(found) != 1:
            sys.exit("Expected one tree for %s, found %d" % (j["id"], len(found)))
        result[j["id"]] = found[0]
    root_job = next(j for j in plan["jobs"] if j["root"])
    tree = parse(open(result[root_job["id"]]).read())
    for j in reversed(plan["jobs"]):
        if j["root"]:
            continue
        sub = parse(open(result[j["id"]]).read())
        part = set(leaves(sub)) - set(j["outgroups"])
        _, cut = find_clade(sub, part, j["outgroups"][0])
        outside = next(t for t in leaves(tree) if t not in j["reps"])
        tree, node = find_clade(tree, j["reps"], outside)
        node.kids, node.name = cut.kids, None
    constraint = parse(open(a.constraint).read())
    broken = check(constraint, tree)
    with open(a.out, "w") as fo:
        fo.write(write(tree, True) + ";\n")
    if broken:
        sys.exit("%d clades of the constraint tree are not monophyletic in the grafted tree" % broken)
    print("grafted %d tips from %d parts" % (len(leaves(tree)), len(plan["jobs"])))


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--version", action="version", version=__version__)
    sub = p.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("split")
    s.add_argument("--alignment", required=True)
    s.add_argument("--constraint", required=True)
    s.add_argument("--partition-tree")
    s.add_argument("--outdir", required=True)
    s.add_argument("--plan", required=True)
    s.add_argument("--max", type=int, default=2000, help="most tips in a part")
    s.add_argument("--min", type=int, default=10, help="fewest tips for which a child is closed as a part")
    s.add_argument("--reps", type=int, default=3, help="tips representing a closed part in its parent")
    s.add_argument("--outgroups", type=int, default=4, help="tips from neighbouring parts added to a part")
    g = sub.add_parser("graft")
    g.add_argument("--plan", required=True)
    g.add_argument("--constraint", required=True)
    g.add_argument("--out", required=True)
    g.add_argument("results", nargs="+")
    a = p.parse_args()
    sys.setrecursionlimit(1_000_000)
    threading.stack_size(512 * 1024 * 1024)
    code = []

    def run():
        try:
            (split if a.cmd == "split" else graft)(a)
        except SystemExit as e:
            code.append(e.code)
        except BaseException:
            traceback.print_exc()
            code.append(1)

    t = threading.Thread(target=run)
    t.start()
    t.join()
    sys.exit(code[0] if code else 0)


if __name__ == "__main__":
    main()
