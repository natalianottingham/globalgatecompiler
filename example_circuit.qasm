// Example circuit used in the tutorial notebook
// 10-Qubit QFT Adder Circuit from Quantum Circuit Benchmarks
// https://github.com/jmbaker94/quantumcircuitbenchmarks

OPENQASM 2.0;
include "qelib1.inc";

qreg q[10];
ry(-pi/2) q[4];
rz(-3*pi/4) q[4];
cx q[4],q[3];
rz(-pi/4) q[3];
cx q[4],q[3];
rz(-3*pi/4) q[3];
ry(pi/2) q[3];
rz(pi/4) q[3];
rz(pi/8) q[4];
cx q[4],q[2];
rz(-pi/8) q[2];
cx q[4],q[2];
rz(pi/8) q[2];
cx q[3],q[2];
rz(-pi/4) q[2];
cx q[3],q[2];
rz(-3*pi/4) q[2];
ry(pi/2) q[2];
rz(pi/4) q[2];
rz(pi/8) q[3];
rz(pi/16) q[4];
cx q[4],q[1];
rz(-pi/16) q[1];
cx q[4],q[1];
rz(pi/16) q[1];
cx q[3],q[1];
rz(-pi/8) q[1];
cx q[3],q[1];
rz(pi/8) q[1];
cx q[2],q[1];
rz(-pi/4) q[1];
cx q[2],q[1];
rz(-3*pi/4) q[1];
ry(pi/2) q[1];
rz(pi/4) q[1];
rz(pi/8) q[2];
rz(pi/16) q[3];
rz(pi/32) q[4];
cx q[4],q[0];
rz(-pi/32) q[0];
cx q[4],q[0];
rz(pi/32) q[0];
cx q[3],q[0];
rz(-pi/16) q[0];
cx q[3],q[0];
rz(pi/16) q[0];
cx q[2],q[0];
rz(-pi/8) q[0];
cx q[2],q[0];
rz(pi/8) q[0];
cx q[1],q[0];
rz(-pi/4) q[0];
cx q[1],q[0];
rz(-3*pi/4) q[0];
ry(pi/2) q[0];
cx q[0],q[4];
cx q[1],q[3];
rz(pi/2) q[2];
cx q[3],q[1];
cx q[1],q[3];
rz(pi/2) q[1];
rz(pi/2) q[3];
cx q[4],q[0];
cx q[0],q[4];
rz(pi/2) q[0];
rz(pi/2) q[4];
cx q[5],q[1];
rz(-pi/2) q[1];
cx q[5],q[1];
rz(0.8862269254527579) q[1];
cx q[5],q[1];
rz(-0.8862269254527579) q[1];
cx q[5],q[1];
rz(0.6656676819001949) q[1];
cx q[5],q[1];
rz(-0.6656676819001949) q[1];
cx q[5],q[1];
rz(0.5769175339249947) q[1];
cx q[5],q[1];
rz(-0.5769175339249947) q[1];
cx q[5],q[1];
rz(0.5370835753981845) q[1];
cx q[5],q[1];
rz(-0.5370835753981845) q[1];
cx q[5],q[1];
cx q[6],q[2];
rz(-pi/2) q[2];
cx q[6],q[2];
rz(0.8862269254527579) q[2];
cx q[6],q[2];
rz(-0.8862269254527579) q[2];
cx q[6],q[2];
rz(0.6656676819001949) q[2];
cx q[6],q[2];
rz(-0.6656676819001949) q[2];
cx q[6],q[2];
rz(0.5769175339249947) q[2];
cx q[6],q[2];
rz(-0.5769175339249947) q[2];
cx q[6],q[2];
rz(0.5370835753981845) q[2];
cx q[6],q[2];
rz(-0.5370835753981845) q[2];
cx q[6],q[2];
rz(-pi/8) q[2];
cx q[7],q[3];
rz(-pi/2) q[3];
cx q[7],q[3];
rz(0.8862269254527579) q[3];
cx q[7],q[3];
rz(-0.8862269254527579) q[3];
cx q[7],q[3];
rz(0.6656676819001949) q[3];
cx q[7],q[3];
rz(-0.6656676819001949) q[3];
cx q[7],q[3];
rz(0.5769175339249947) q[3];
cx q[7],q[3];
rz(-0.5769175339249947) q[3];
cx q[7],q[3];
rz(0.5370835753981845) q[3];
cx q[7],q[3];
rz(-0.5370835753981845) q[3];
cx q[7],q[3];
cx q[1],q[3];
cx q[3],q[1];
cx q[1],q[3];
rz(-pi/4) q[1];
rz(-pi/16) q[3];
cx q[8],q[4];
rz(-pi/2) q[4];
cx q[8],q[4];
rz(0.8862269254527579) q[4];
cx q[8],q[4];
rz(-0.8862269254527579) q[4];
cx q[8],q[4];
rz(0.6656676819001949) q[4];
cx q[8],q[4];
rz(-0.6656676819001949) q[4];
cx q[8],q[4];
rz(0.5769175339249947) q[4];
cx q[8],q[4];
rz(-0.5769175339249947) q[4];
cx q[8],q[4];
rz(0.5370835753981845) q[4];
cx q[8],q[4];
rz(-0.5370835753981845) q[4];
cx q[8],q[4];
cx q[9],q[0];
rz(-pi/2) q[0];
cx q[9],q[0];
rz(0.8862269254527579) q[0];
cx q[9],q[0];
rz(-0.8862269254527579) q[0];
cx q[9],q[0];
rz(0.6656676819001949) q[0];
cx q[9],q[0];
rz(-0.6656676819001949) q[0];
cx q[9],q[0];
rz(0.5769175339249947) q[0];
cx q[9],q[0];
rz(-0.5769175339249947) q[0];
cx q[9],q[0];
rz(0.5370835753981845) q[0];
cx q[9],q[0];
rz(-0.5370835753981845) q[0];
cx q[9],q[0];
cx q[0],q[4];
cx q[4],q[0];
cx q[0],q[4];
ry(pi/2) q[0];
rx(pi) q[0];
cx q[1],q[0];
rz(pi/4) q[0];
cx q[1],q[0];
rz(-pi/4) q[0];
ry(pi/2) q[1];
rx(pi) q[1];
cx q[2],q[0];
rz(pi/8) q[0];
cx q[2],q[0];
rz(-pi/8) q[0];
rz(-pi/4) q[2];
cx q[2],q[1];
rz(pi/4) q[1];
cx q[2],q[1];
rz(-pi/4) q[1];
ry(pi/2) q[2];
rx(pi) q[2];
cx q[3],q[0];
rz(pi/16) q[0];
cx q[3],q[0];
rz(-pi/16) q[0];
rz(-pi/8) q[3];
cx q[3],q[1];
rz(pi/8) q[1];
cx q[3],q[1];
rz(-pi/8) q[1];
rz(-pi/4) q[3];
cx q[3],q[2];
rz(pi/4) q[2];
cx q[3],q[2];
rz(-pi/4) q[2];
ry(pi/2) q[3];
rx(pi) q[3];
rz(-pi/32) q[4];
cx q[4],q[0];
rz(pi/32) q[0];
cx q[4],q[0];
rz(-pi/32) q[0];
rz(-pi/16) q[4];
cx q[4],q[1];
rz(pi/16) q[1];
cx q[4],q[1];
rz(-pi/16) q[1];
rz(-pi/8) q[4];
cx q[4],q[2];
rz(pi/8) q[2];
cx q[4],q[2];
rz(-pi/8) q[2];
rz(-pi/4) q[4];
cx q[4],q[3];
rz(pi/4) q[3];
cx q[4],q[3];
rz(-pi/4) q[3];
ry(pi/2) q[4];
rx(pi) q[4];