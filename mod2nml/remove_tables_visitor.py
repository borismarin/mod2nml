import nmodl
from nmodl import to_nmodl as N

class RemoveTablesVisitor(nmodl.visitor.AstVisitor):
    def visit_statement_block(self, node):
        stmts = node.statements
        for i,s in enumerate(node.statements):
            if s.is_table_statement():  # can't use visit_table statement
                s_as_string = nmodl.dsl.ast.String('; ' + nmodl.to_nmodl(s))
                stmts[i] = nmodl.dsl.ast.LineComment(s_as_string)
        node.statements = stmts

