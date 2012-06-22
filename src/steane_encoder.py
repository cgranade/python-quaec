import qecc as q
import operator as op
from itertools import starmap, product
from collections import defaultdict

def error_to_pauli(error):
    if error == q.I.as_clifford():
        return "I"
    if error == q.X.as_clifford():
        return "X"
    if error == q.Y.as_clifford():
        return "Y"
    if error == q.Z.as_clifford():
        return "Z"

if __name__ == "__main__":
    stab_code = q.StabilizerCode.steane_code()
    enc = stab_code.encoding_cliffords().next()
    dec = enc.inv()

    print "Encoder circuit:"
    print "========================================================================"
    circuit = enc.as_bsm().circuit_decomposition().replace_cz_by_cnot().cancel_selfinv_gates()
    circuit = sum(circuit.group_by_time(pad_with_waits=True), q.Circuit())

    # We use group_by_time to find where to add wait locations, then
    # recombine the circuit to obtain the decoder.
    print circuit


    print "\n\n"
    print "Syndrome propagation:"
    print "========================================================================"

    faults = [q.elem_gen(7, idx, P) for idx in range(7) for P in ['X', 'Y', 'Z']]
    synd_meas = [q.elem_gen(7, idx, kind) for idx, kind in zip(range(1,7), 'ZZZZZZ')]

    stab = q.StabilizerCode.unencoded_state(nq_logical=1, nq_ancilla=6)

    print "Initial stabilizer code:"
    print stab

    print "Encoding:"
    print enc(stab)


    print "Using syndrome measurement operators:"
    for idx, meas in enumerate(synd_meas):
        print '\t{}\t{}'.format(idx, meas.str_sparse(incl_ph=False))
    print ''
        
    recovery = defaultdict(lambda: q.I.as_clifford())
        
    for fault in faults:
        eff = dec * fault.as_clifford() * enc
        syndrome = tuple([eff.conjugate_pauli(meas).ph / 2 for meas in synd_meas])
        error = q.Clifford([eff.xout[0][0]], [eff.zout[0][0]])
        if syndrome not in recovery:
            recovery[syndrome] = error
        else:
            if recovery[syndrome] != error:
                raise RuntimeError("Collision.")

    print "Recovery operators for all weight-1 errors:"

    inv_recovery = defaultdict(list)
    for syndrome in product(range(2), repeat=6):
        inv_recovery[error_to_pauli(recovery[syndrome])].append(syndrome)
        
    for recovery_operator in ['I', 'X', 'Y', 'Z']:
        syndromes = inv_recovery[recovery_operator]
        print "{}: {}".format(recovery_operator, ", ".join([
            "".join(map(str, syndrome)) for syndrome in syndromes
        ]))
        
