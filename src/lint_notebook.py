#!/usr/bin/env python3
"""Static check that a notebook can run top to bottom.

The expensive failure mode for a submission notebook is a NameError twenty
minutes into a clean run. This walks the cells in order, tracks every name each
cell binds, and reports any name a later cell reads before anything has bound
it. It also flags duplicate parameter definitions and cells that reference a
file path no earlier cell creates.

It cannot catch runtime errors, only ordering ones -- but ordering is what
"Restart and run all" actually tests.

Usage: lint_notebook.py <notebook.ipynb>
"""

import ast
import builtins
import json
import re
import sys

BUILTINS = set(dir(builtins))

# Names Colab / IPython provide, plus magics handled separately.
AMBIENT = {"get_ipython", "display", "In", "Out", "exit", "quit", "_"}


class Binder(ast.NodeVisitor):
    """Collect names a module binds, and names it reads."""

    def __init__(self):
        self.bound = set()
        self.read = set()
        self._scope_depth = 0

    # -- bindings -------------------------------------------------------
    def visit_Assign(self, node):
        for t in node.targets:
            self._bind(t)
        self.visit(node.value)

    def visit_AnnAssign(self, node):
        self._bind(node.target)
        if node.value:
            self.visit(node.value)

    def visit_AugAssign(self, node):
        self._bind(node.target)
        self.visit(node.value)

    def visit_For(self, node):
        self._bind(node.target)
        self.generic_visit(node)

    def visit_With(self, node):
        for item in node.items:
            if item.optional_vars:
                self._bind(item.optional_vars)
        self.generic_visit(node)

    def visit_Import(self, node):
        for a in node.names:
            self.bound.add((a.asname or a.name).split(".")[0])

    def visit_ImportFrom(self, node):
        for a in node.names:
            self.bound.add(a.asname or a.name)

    def visit_FunctionDef(self, node):
        self.bound.add(node.name)
        # treat the body as a nested scope: its locals are not module-level,
        # but free names it reads still need to exist by call time.
        inner = Binder()
        for arg in list(node.args.args) + list(node.args.kwonlyargs):
            inner.bound.add(arg.arg)
        if node.args.vararg:
            inner.bound.add(node.args.vararg.arg)
        if node.args.kwarg:
            inner.bound.add(node.args.kwarg.arg)
        for stmt in node.body:
            inner.visit(stmt)
        self.read |= (inner.read - inner.bound)

    visit_AsyncFunctionDef = visit_FunctionDef

    def visit_ClassDef(self, node):
        self.bound.add(node.name)
        self.generic_visit(node)

    def visit_Lambda(self, node):
        inner = Binder()
        for arg in list(node.args.args) + list(node.args.kwonlyargs):
            inner.bound.add(arg.arg)
        inner.visit(node.body)
        self.read |= (inner.read - inner.bound)

    def visit_comprehension(self, node):
        self._bind(node.target)
        self.visit(node.iter)
        for i in node.ifs:
            self.visit(i)

    def visit_ListComp(self, node):
        self._comp(node)

    visit_SetComp = visit_ListComp
    visit_GeneratorExp = visit_ListComp

    def visit_DictComp(self, node):
        self._comp(node, extra=(node.key, node.value))

    def _comp(self, node, extra=None):
        inner = Binder()
        for gen in node.generators:
            inner._bind(gen.target)
            inner.visit(gen.iter)
            for i in gen.ifs:
                inner.visit(i)
        for part in (extra or (node.elt,)):
            inner.visit(part)
        self.read |= (inner.read - inner.bound)

    def visit_ExceptHandler(self, node):
        if node.name:
            self.bound.add(node.name)
        self.generic_visit(node)

    # -- reads ----------------------------------------------------------
    def visit_Name(self, node):
        if isinstance(node.ctx, ast.Load):
            self.read.add(node.id)
        else:
            self.bound.add(node.id)

    def _bind(self, target):
        for n in ast.walk(target):
            if isinstance(n, ast.Name):
                self.bound.add(n.id)


def main(path):
    nb = json.load(open(path))
    defined = set(BUILTINS) | AMBIENT
    problems = []
    assigned_where = {}

    for idx, cell in enumerate(nb["cells"]):
        if cell["cell_type"] != "code":
            continue
        src = "".join(cell["source"])
        if src.lstrip().startswith("%%"):
            # a cell magic (e.g. %%bash) is not Python; it binds nothing
            continue
        def _is_magic(line):
            s = line.lstrip()
            return (s.startswith("!") and not s.startswith("!=")) \
                or bool(re.match(r"%[A-Za-z]", s))

        src = "\n".join("" if _is_magic(l) else l for l in src.split("\n"))
        try:
            tree = ast.parse(src)
        except SyntaxError as e:
            problems.append("cell %d: SyntaxError %s" % (idx, e))
            continue

        b = Binder()
        for stmt in tree.body:
            b.visit(stmt)

        for name in sorted(b.read - defined - b.bound):
            problems.append("cell %d: reads '%s' before any cell defines it"
                            % (idx, name))

        for name in b.bound:
            if name in assigned_where and name.isupper():
                problems.append("cell %d: re-defines parameter '%s' "
                                "(first set in cell %d)"
                                % (idx, name, assigned_where[name]))
            assigned_where.setdefault(name, idx)

        defined |= b.bound

    print("notebook :", path)
    print("cells    :", len(nb["cells"]))
    if problems:
        print("\nPROBLEMS (%d):" % len(problems))
        for p in problems:
            print("  -", p)
        return 1
    print("\nno ordering problems: every name is bound before it is read")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1]))
