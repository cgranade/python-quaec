import qecc as q
import operator as op

stab_code = q.StabilizerCode.steane_code()
enc = stab_code.encoding_cliffords().next()
dec = enc.as_bsm().inv().as_clifford()


circuit = enc.as_bsm().circuit_decomposition()
print '======='
circuit.replace_cz_by_cnot()
print circuit
print '======='
circuit.cancel_selfinv_gates()
print circuit
print '======='

# We use group_by_time to find where to add wait locations, then
# recombine the circuit to obtain the decoder.
print sum(circuit.group_by_time(pad_with_waits=True), q.Circuit())
print len(circuit), circuit.depth

faults = [q.elem_gen(7, idx, P) for idx in range(7) for P in ['X', 'Z']]
for fault in faults:
    print "Fault {}:".format(fault)
    #print fault.as_clifford()
    eff = enc * fault.as_clifford() * dec
    print "Syndrome: [{}]".format(
        ", ".join(
            [str(0 if eff.xout[idx].ph == 0 else 1) for idx in range(1,4)] +
            [str(0 if eff.zout[idx].ph == 0 else 1) for idx in range(4,7)]
        )
    )
    print ''
    print "Effective single-qubit Clifford:"
    print q.Clifford([eff.xout[0][0]], [eff.zout[0][0]])
    print ''
    print '--------'
    
