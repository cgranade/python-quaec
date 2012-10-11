import qecc as q

def test_bsm_reductions():
    cnot_in = map(q.Pauli, [
        'XX', 'IX', 'ZI', 'ZZ'
    ])

    cnot_out = map(q.Pauli, [
        'XI', 'IX', 'ZI', 'IZ'
    ])

    C = q.gen_cliff(cnot_in, cnot_out)
    testbsm=q.cnot(2,0,1).as_bsm()
    print testbsm
    a=testbsm.left_H(0)
    print a
    b=testbsm.left_SWAP(0,1)
    print b

    return b == (q.swap(2, 0, 1) * q.hadamard(2, 0) * C).as_bsm()

if __name__ == "__main__":
    assert test_bsm_reductions()
