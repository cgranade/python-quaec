import qecc as q
import operator as op
from itertools import starmap

stab_code = q.StabilizerCode.steane_code()
enc = stab_code.encoding_cliffords().next()
dec = enc.as_bsm().inv().as_clifford()

print "Encoder circuit:"
print "========================================================================"
circuit = enc.as_bsm().circuit_decomposition()
print '======='
circuit.replace_cz_by_cnot()
print circuit
print '======='
circuit.cancel_selfinv_gates()
print circuit

# We use group_by_time to find where to add wait locations, then
# recombine the circuit to obtain the decoder.
print sum(circuit.group_by_time(pad_with_waits=True), q.Circuit())
print len(circuit), circuit.depth


print "\n\n"
print "Syndrome propagation:"
print "========================================================================"

faults = [q.elem_gen(7, idx, P) for idx in range(7) for P in ['X', 'Z']]
synd_meas = [q.elem_gen(7, idx, kind) for idx, kind in zip(range(1,7), 'XXXZZZ')]

print "Using syndrome measurement operators:"
for idx, meas in enumerate(synd_meas):
    print '\t{}\t{}'.format(idx, meas.str_sparse(incl_ph=False))
print ''
    
for fault in faults:
    eff = enc * fault.as_clifford() * dec
    print "Fault {} yields syndrome [{}].".format(
        fault.str_sparse(incl_ph=False),
        ", ".join(str(eff.conjugate_pauli(meas).ph / 2) for meas in synd_meas)
        )
    print ''
    print "Effective single-qubit Clifford:"
    print q.Clifford([eff.xout[0][0]], [eff.zout[0][0]])
    print '--------'
    
