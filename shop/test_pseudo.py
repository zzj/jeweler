import pseudo
from nose.tools import *

def test_pseudo_set():
    ps = pseudo.PseudoSet()
    assert_equal(True, ps.is_pseudo_gene("Gm10557"))
    assert_equal(False, ps.is_pseudo_gene("Noname"))
