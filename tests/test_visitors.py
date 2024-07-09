import nmodl.dsl as nmodl
from mod2nml.rename_visitor import RenameReusedVisitor
from mod2nml.remove_tables_visitor import RemoveTablesVisitor


def test_rename_reused():
    mod = """
    DERIVATIVE states   {
        a = 10
        b = 20
        z = 2 * a * f(a)
        a = 30
        cai' = a + b + g(a^2)
        a = 1
        cai' = a + h(a^a)
    }
    """
    driver = nmodl.NmodlDriver()
    modast = driver.parse_string(mod)

    param_visitor = RenameReusedVisitor()
    modast.accept(param_visitor)

    print(nmodl.to_nmodl(modast))


def test_remove_unsopported():
    mod = """
    PROCEDURE settables(v) {
      TABLE minf, mtau FROM -120 TO 40 WITH 641

      minf  = 1 / ( 1 + exp( ( - v - 48 ) / 10 ) )

      if( v < -40.0 ) {
              mtau = 0.025 + 0.14 * exp( ( v + 40 ) / 10 )
      }else{
              mtau = 0.02 + 0.145 * exp( ( - v - 40 ) / 10 )
      }
    }
    """
    driver = nmodl.NmodlDriver()
    modast = driver.parse_string(mod)

    procstmts = modast.blocks[0].statement_block.statements
    print(type(procstmts[0]), f'"{procstmts[0]}"')
    assert procstmts[0].is_table_statement

    tables_visitor = RemoveTablesVisitor()
    modast.accept(tables_visitor)

    procstmts = modast.blocks[0].statement_block.statements
    print(type(procstmts[0]), f'"{procstmts[0]}"')
    assert procstmts[0].is_line_comment

    #print(nmodl.to_nmodl(modast))

