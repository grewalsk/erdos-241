"""Standalone B_3 verifier for Erdős Problem #241.

Convention (fixed for this whole project):
  A finite set A of integers is a B_3 set if every integer has AT MOST ONE
  representation as a+b+c with a <= b <= c, a,b,c in A.  Repetition is
  allowed (a+a+b and a+a+a are representations); two representations are the
  same iff they are equal as multisets.

This module is deliberately brute-force and independent of every solver /
construction in this project: it enumerates ALL multisets of size 3 via
itertools.combinations_with_replacement and checks the sums for collisions.
It is the ground truth everything else is checked against.

Also provides the modular variant (sums taken mod m) and a B_2 (Sidon)
checker used for sanity tests.

Usage:
    python3 verify.py 1 2 5 11 22 33          # check a set, integer B_3
    python3 verify.py --mod 63 0 1 5 21       # check a set, B_3 in Z_63
"""

from itertools import combinations_with_replacement
import sys


def b3_violations(A, modulus=None, max_report=20):
    """Return a list of collisions ((multiset1, multiset2, sum)).

    Empty list  <=>  A is a B_3 set (in Z if modulus is None, else in Z_modulus).
    Raises ValueError if A contains duplicate elements (it must be a SET;
    in the modular case, distinct residues mod modulus).
    """
    if modulus is None:
        elems = sorted(A)
    else:
        elems = sorted(a % modulus for a in A)
    if len(elems) != len(set(elems)):
        raise ValueError("input contains duplicate elements (after reduction)")
    seen = {}
    viol = []
    for t in combinations_with_replacement(elems, 3):
        s = sum(t)
        if modulus is not None:
            s %= modulus
        if s in seen:
            viol.append((seen[s], t, s))
            if len(viol) >= max_report:
                return viol
        else:
            seen[s] = t
    return viol


def is_b3(A, modulus=None):
    """True iff A is a B_3 set (multiset convention)."""
    return not b3_violations(A, modulus=modulus, max_report=1)


def b2_violations(A, modulus=None, max_report=20):
    """Sidon (B_2) violations, same conventions, for sanity checks."""
    if modulus is None:
        elems = sorted(A)
    else:
        elems = sorted(a % modulus for a in A)
    if len(elems) != len(set(elems)):
        raise ValueError("input contains duplicate elements (after reduction)")
    seen = {}
    viol = []
    for t in combinations_with_replacement(elems, 2):
        s = sum(t)
        if modulus is not None:
            s %= modulus
        if s in seen:
            viol.append((seen[s], t, s))
            if len(viol) >= max_report:
                return viol
        else:
            seen[s] = t
    return viol


def is_b2(A, modulus=None):
    return not b2_violations(A, modulus=modulus, max_report=1)


def _selftest():
    # Known tiny facts, hand-checkable.
    assert is_b3([1])
    assert is_b3([1, 2])            # sums 3,4,5,6 distinct
    assert not is_b3([1, 2, 3])     # 1+2+3 = 6 = 2+2+2
    assert is_b3([1, 2, 5])         # sums 3,4,7,5,8,11,6,9,12,15 distinct
    # {1,2,5,11} is NOT B_3: 2+2+11 = 15 = 5+5+5
    assert not is_b3([1, 2, 5, 11])
    # {1,2,5,14} is B_3 (hand-checked: new sums 16,17,18,20,21,24,29,30,33,42)
    assert is_b3([1, 2, 5, 14])
    # B_3 => B_2 (padding argument); spot-check on a few sets
    for S in ([1, 2, 5, 14], [1, 2, 5], [3, 4, 10, 21]):
        if is_b3(S):
            assert is_b2(S), f"B_3 set {S} not Sidon?!"
    # modular: {0,1,3} in Z_7?  triple sums mod 7: 0,1,3,2,4,6,3->collision(0+0+3 vs 1+1+1)
    assert not is_b3([0, 1, 3], modulus=7)
    assert is_b3([0, 1], modulus=7)
    # a set that is B_3 over Z can fail mod small m
    assert is_b3([1, 2, 5])
    assert not is_b3([1, 2, 5], modulus=9)  # 1+1+5=7, 2+2+5=9=0, 5+5+5=15=6, 1+2+2=5, ... check: brute decides
    print("verify.py selftest OK")


if __name__ == "__main__":
    args = sys.argv[1:]
    if not args:
        _selftest()
        sys.exit(0)
    modulus = None
    if args[0] == "--mod":
        modulus = int(args[1])
        args = args[2:]
    A = [int(x) for x in args]
    v = b3_violations(A, modulus=modulus)
    if not v:
        print(f"OK: B_3 set of size {len(A)}" + (f" in Z_{modulus}" if modulus else ""))
    else:
        print(f"FAIL: {len(v)} collision(s) shown:")
        for t1, t2, s in v:
            print(f"  {t1} and {t2} both sum to {s}")
        sys.exit(1)
