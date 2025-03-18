import { Bindings } from 'o1js';

import {
  Point,
  TwistedCurve,
  simpleMapToCurve,
  Field3,
} from './provable/twisted-curve.js';

export { TwistedCurve, Point, simpleMapToCurve, Field3 };

import {
  AffineTwistedCurve,
  GroupAffineTwisted,
  affineTwistedAdd,
  affineTwistedDouble,
  affineTwistedZero,
  TwistedCurveParams,
} from './crypto/elliptic-curve.js';

export {
  AffineTwistedCurve,
  GroupAffineTwisted,
  affineTwistedAdd,
  affineTwistedDouble,
  affineTwistedZero,
  TwistedCurveParams,
};

// Parameters used in Ed25519 (EdDSA algorithm for edwards25519 curve)
// https://datatracker.ietf.org/doc/html/rfc8032#section-5.1
const edwards25519Params: TwistedCurveParams = {
  name: 'edwards25519',
  modulus: Bindings.exampleFields.f25519.modulus, // 2^255 - 19
  order: 0x1000000000000000000000000000000014def9dea2f79cd65812631a5cf5d3edn, //2^252 + 27742317777372353535851937790883648493,
  cofactor: 8n,
  generator: {
    x: 0x216936d3cd6e53fec0a4e231fdd6dc5c692cc7609525a7b2c9562d608f25d51an, // <=> 15112221349535400772501151409588531511454012693041857206046113283949847762202
    y: 0x6666666666666666666666666666666666666666666666666666666666666658n, // <=> 4/5 mod p <=> 46316835694926478169428394003475163141307993866256225615783033603165251855960
  },
  a: 0x7fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffecn, // <=> -1 mod p <=> 57896044618658097711785492504343953926634992332820282019728792003956564819948
  d: 0x52036cee2b6ffe738cc740797779e89800700a4d4141d8ab75eb4dca135978a3n, // -121665/121666 mod p <=> 37095705934669439343138083508754565189542113879843219016388785533085940283555
};

const TwistedCurveParams = {
  Edwards25519: edwards25519Params,
};
