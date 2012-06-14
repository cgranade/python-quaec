import qecc as q
import operator as op

stab_in = q.PauliList([q.elem_gen(7, idx, 'Z') for idx in range(1,7)])
stab_out = q.PauliList('ZZZZIII', 'ZZIIZZI', 'ZIZIZIZ', 'XXXXIII', 'XXIIXXI', 'XIXIXIX')

Xbars_in = q.PauliList(['XIIIIII'] + [q.Unspecified] * 6)
Zbars_in = q.PauliList(['ZIIIIII'] + stab_in)

Xbars_out = q.PauliList(['XXXXXXX'] + [q.Unspecified] * 6)
Zbars_out = q.PauliList(['ZZZZZZZ'] + stab_out)

C_in = q.Clifford(Xbars_in, Zbars_in).constraint_completions().next()
C_out = q.Clifford(Xbars_out, Zbars_out).constraint_completions().next()
enc = C_in.as_bsm().inv().as_clifford() * C_out

circuit = enc.as_bsm().circuit_decomposition()
print type(circuit)
assert isinstance(circuit, q.Circuit)
print circuit
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
print len(circuit), circuit.n_timesteps
