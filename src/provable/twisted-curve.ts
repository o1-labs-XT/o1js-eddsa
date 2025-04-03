import {
  Provable,
  assert,
  Gadgets,
  provable,
  Core,
  Field,
  UInt8,
  Bool,
  AlmostForeignField,
  Bytes,
  ProvablePureExtended,
  createForeignField,
  provableFromClass,
} from 'o1js';

import {
  AffineTwistedCurve,
  GroupAffineTwisted,
  affineTwistedAdd,
  affineTwistedDouble,
  affineTwistedZero,
  createAffineTwistedCurve,
} from '../crypto/elliptic-curve.js';
import { TwistedCurveParams } from '../index.js';

export {
  TwistedCurve,
  Point,
  simpleMapToCurve,
  Field3,
  scale,
  Eddsa,
  toPoint,
  FlexiblePoint,
  encode,
  ForeignTwisted,
  createForeignTwisted,
  TwistedCurves,
};

const { Field3, arrayGetGeneric, ForeignField, SHA2, multiRangeCheck } =
  Gadgets;

const { l2Mask } = Gadgets.Constants;
const { mod } = Core.FiniteField;

const TwistedCurve = {
  add,
  double,
  negate,
  assertOnCurve,
  scale,
  multiScalarMul,
  assertInSubgroup,
};

type Field3 = [Field, Field, Field];

type FlexiblePoint = {
  x: AlmostForeignField | Field3 | bigint | number;
  y: AlmostForeignField | Field3 | bigint | number;
};

function toPoint({ x, y }: ForeignTwisted): Point {
  return { x: x.value, y: y.value };
}

/**
 * Non-zero twisted elliptic curve point.
 */
type Point = { x: Field3; y: Field3 };
type point = { x: bigint; y: bigint };

const Point = {
  from({ x, y }: point): Point {
    return { x: Field3.from(x), y: Field3.from(y) };
  },
  toBigint({ x, y }: Point) {
    let x_ = Field3.toBigint(x);
    let y_ = Field3.toBigint(y);
    return { x: x_, y: y_ };
  },
  isConstant: (P: Point) => Provable.isConstant(Point, P),

  /**
   * Random point on the curve.
   */
  random(Curve: AffineTwistedCurve) {
    return Point.from(random(Curve));
  },

  provable: provable({ x: Field3, y: Field3 }),
  /**
   * On input a compressed representation of a Edwards25519 point as 32 bytes of
   * hexadecimal string, return the point with bigint coordinates.
   *
   * @param hex 32 bytes of hexadecimal string
   * @returns point with bigint coordinates {x, y}
   */
  fromHex(hex: string): point {
    const y = BigInt(`0x${hex}`) & ((BigInt(1) << 255n) - 1n); // y (mask top bit)
    const x_0 = (BigInt(`0x${hex}`) >> 255n) & 1n; // parity bit for x

    if (y >= Curve.modulus) {
      throw new Error(`Invalid y value: ${y} is larger tan the field size.`);
    }

    let x = recoverX(y, x_0);

    return { x, y };
  },
};

function add(
  p1: Point,
  p2: Point,
  Curve: { modulus: bigint; a: bigint; d: bigint }
) {
  let { x: x1, y: y1 } = p1;
  let { x: x2, y: y2 } = p2;
  let f = Curve.modulus;
  let a = Curve.a;
  let d = Curve.d;

  // constant case
  if (Point.isConstant(p1) && Point.isConstant(p2)) {
    let p3 = affineTwistedAdd(Point.toBigint(p1), Point.toBigint(p2), f, a, d);
    return Point.from(p3);
  }

  assert(
    Curve.modulus > l2Mask + 1n,
    'Base field moduli smaller than 2^176 are not supported'
  );

  // the formula for point addition is well defined for curves in use,
  // so we don't need to check that the denominators are non-zero

  // x3 = (x1 * y2 + y1 * x2) / (1 + d * x1 * x2 * y1 * y2)
  // y3 = (y1 * y2 - a * x1 * x2) / (1 - d * x1 * x2 * y1 * y2)

  let x1x2 = ForeignField.mul(x1, x2, f);
  let y1y2 = ForeignField.mul(y1, y2, f);
  let x1y2 = ForeignField.mul(x1, y2, f);
  let y1x2 = ForeignField.mul(y1, x2, f);
  let ax1x2 = ForeignField.mul(Field3.from(a), x1x2, f);

  let x3Num = ForeignField.add(x1y2, y1x2, f);
  let y3Num = ForeignField.sub(y1y2, ax1x2, f);

  let x1x2y1y2 = ForeignField.mul(x1x2, y1y2, f);
  let dx1x2y1y2 = ForeignField.mul(Field3.from(d), x1x2y1y2, f);

  let one = Field3.from(1n);
  let x3Denom = ForeignField.add(one, dx1x2y1y2, f);
  let y3Denom = ForeignField.sub(one, dx1x2y1y2, f);

  let x3 = ForeignField.div(x3Num, x3Denom, f);
  let y3 = ForeignField.div(y3Num, y3Denom, f);

  ForeignField.assertAlmostReduced(
    [x1x2, y1y2, x3Num, y3Num, x1x2y1y2, x3Denom, y3Denom, x3, y3],
    f
  );

  return { x: x3, y: y3 };
}

function double(
  p1: Point,
  Curve: { modulus: bigint; a: bigint; d: bigint }
): Point {
  let { x: x1, y: y1 } = p1;
  let f = Curve.modulus;
  let d = Curve.d;

  // constant case
  if (Point.isConstant(p1)) {
    let p3 = affineTwistedDouble(Point.toBigint(p1), f, Curve.a, Curve.d);
    return Point.from(p3);
  }

  // x3 = 2*x1*y1 / (1 + d * x1^2 * y1^2)
  // y3 = (y1^2 - a * x1^2) / (1 - d * x1^2 * y1^2)
  let one = Field3.from(1n);
  let a = Field3.from(Curve.a);
  let x1x1 = ForeignField.mul(x1, x1, f);
  let y1y1 = ForeignField.mul(y1, y1, f);
  let x1y1 = ForeignField.mul(x1, y1, f);
  let ax1x1 = ForeignField.mul(a, x1x1, f);
  let x3Num = ForeignField.add(x1y1, x1y1, f);
  let y3Num = ForeignField.sub(y1y1, ax1x1, f);
  let x1x1y1y1 = ForeignField.mul(x1x1, y1y1, f);
  let dx1x1y1y1 = ForeignField.mul(Field3.from(d), x1x1y1y1, f);
  let x3Den = ForeignField.add(one, dx1x1y1y1, f);
  let y3Den = ForeignField.sub(one, dx1x1y1y1, f);
  let x3 = ForeignField.div(x3Num, x3Den, f);
  let y3 = ForeignField.div(y3Num, y3Den, f);

  ForeignField.assertAlmostReduced([x3Num, y3Num, x3Den, y3Den, x3, y3], f);

  return { x: x3, y: y3 };
}

function negate({ x, y }: Point, Curve: { modulus: bigint }) {
  return { x: ForeignField.neg(x, Curve.modulus), y };
}

function assertOnCurve(
  p: Point,
  { modulus: f, a, d }: { modulus: bigint; a: bigint; d: bigint }
) {
  let { x, y } = p;
  let one = Field3.from(1n);

  // a * x^2 + y^2 = 1 + d * x^2 * y^2

  let x2 = ForeignField.mul(x, x, f);
  let y2 = ForeignField.mul(y, y, f);

  let aTimesX2PlusY2 = ForeignField.add(
    ForeignField.mul(Field3.from(a), x2, f),
    y2,
    f
  );

  let aTimesX2PlusY2Minus1 = ForeignField.sub(aTimesX2PlusY2, one, f);
  let dTimesX2 = ForeignField.mul(Field3.from(d), x2, f);

  ForeignField.assertAlmostReduced([x2, x, y], f);
  ForeignField.assertAlmostReduced([y2, aTimesX2PlusY2Minus1, dTimesX2], f);

  let message: string | undefined;
  if (Point.isConstant(p)) {
    message = `assertOnCurve(): (${x}, ${y}) is not on the curve.`;
  }
  ForeignField.assertMul(dTimesX2, y2, aTimesX2PlusY2Minus1, f, message);
}

/**
 * Twisted curve scalar multiplication, `scalar*point`
 */
function scale(
  scalar: Field3,
  point: Point,
  Curve: AffineTwistedCurve,
  config?: {
    mode?: 'assert-zero' | 'assert-nonzero';
    windowSize?: number;
    multiples?: Point[];
  }
) {
  config = config ?? {};
  config.windowSize ??= Point.isConstant(point) ? 4 : 3;
  return multiScalarMul([scalar], [point], Curve, [config], config.mode);
}

// check whether a point equals a constant point
function equals(p1: Point, p2: point, Curve: { modulus: bigint }) {
  let xEquals = ForeignField.equals(p1.x, p2.x, Curve.modulus);
  let yEquals = ForeignField.equals(p1.y, p2.y, Curve.modulus);
  return xEquals.and(yEquals);
}

// checks whether the twisted elliptic curve point g is in the subgroup defined by [order]g = 0
function assertInSubgroup(g: Point, Curve: AffineTwistedCurve) {
  if (!Curve.hasCofactor) return;
  scale(Field3.from(Curve.order), g, Curve, { mode: 'assert-zero' });
}

function multiScalarMulConstant(
  scalars: Field3[],
  points: Point[],
  Curve: AffineTwistedCurve
): Point {
  let n = points.length;
  assert(scalars.length === n, 'Points and scalars lengths must match');
  assertPositiveInteger(n, 'Expected at least 1 point and scalar');

  // TODO dedicated MSM
  let s = scalars.map(Field3.toBigint);
  let P = points.map(Point.toBigint);
  let sum: GroupAffineTwisted = Curve.zero;
  for (let i = 0; i < n; i++) {
    sum = Curve.add(sum, Curve.scale(P[i], s[i]));
  }
  return Point.from(sum);
}

/**
 * Multi-scalar multiplication:
 *
 * s_0 * P_0 + ... + s_(n-1) * P_(n-1)
 *
 * where P_i are any points.
 *
 * Implementation: We double all points together and leverage a precomputed table of size 2^c to avoid all but every cth addition.
 *
 * Note: this algorithm targets a small number of points
 *
 * TODO: could use lookups for picking precomputed multiples, instead of O(2^c) provable switch
 */
function multiScalarMul(
  scalars: Field3[],
  points: Point[],
  Curve: AffineTwistedCurve,
  tableConfigs: (
    | { windowSize?: number; multiples?: Point[] }
    | undefined
  )[] = [],
  mode?: 'assert-zero' | 'assert-nonzero'
): Point {
  let n = points.length;
  assert(scalars.length === n, 'Points and scalars lengths must match');
  assertPositiveInteger(n, 'Expected at least 1 point and scalar');

  // constant case
  if (scalars.every(Field3.isConstant) && points.every(Point.isConstant)) {
    return multiScalarMulConstant(scalars, points, Curve);
  }

  // parse or build point tables
  let windowSizes = points.map((_, i) => tableConfigs[i]?.windowSize ?? 1);
  let tables = points.map((P, i) =>
    getPointTable(Curve, P, windowSizes[i], tableConfigs[i]?.multiples)
  );

  let maxBits = Curve.Scalar.sizeInBits;

  // slice scalars
  let scalarChunks = scalars.map((s, i) =>
    ForeignField.sliceField3(s, { maxBits, chunkSize: windowSizes[i] })
  );

  // soundness follows because add() and double() are sound, on all inputs that
  // are valid non-zero curve points
  let sum = Point.from(Curve.zero);

  for (let i = maxBits - 1; i >= 0; i--) {
    // add in multiple of each point
    for (let j = 0; j < n; j++) {
      let windowSize = windowSizes[j];
      if (i % windowSize === 0) {
        // pick point to add based on the scalar chunk
        let sj = scalarChunks[j][i / windowSize];
        let sjP =
          windowSize === 1
            ? points[j]
            : arrayGetGeneric(Point.provable, tables[j], sj);

        // ec addition
        sum = add(sum, sjP, Curve);
      }
    }

    if (i === 0) break;

    // jointly double all points
    // (note: the highest couple of bits will not create any constraints because
    // sum is constant; no need to handle that explicitly)
    sum = double(sum, Curve);
  }

  let isZero = equals(sum, affineTwistedZero, Curve);
  if (mode == 'assert-nonzero') {
    isZero.assertFalse();
  } else if (mode == 'assert-zero') {
    isZero.assertTrue();
  }

  return sum;
}

/**
 * Given a point P, create the list of multiples [0, P, 2P, 3P, ..., (2^windowSize-1) * P].
 * This method is provable, but won't create any constraints given a constant point.
 */
function getPointTable(
  Curve: AffineTwistedCurve,
  P: Point,
  windowSize: number,
  table?: Point[]
): Point[] {
  assertPositiveInteger(windowSize, 'invalid window size');
  let n = 1 << windowSize; // n >= 2

  assert(table === undefined || table.length === n, 'invalid table');
  if (table !== undefined) return table;

  table = [Point.from(Curve.zero), P];
  if (n === 2) return table;

  let Pi = double(P, Curve);
  table.push(Pi);
  for (let i = 3; i < n; i++) {
    Pi = add(Pi, P, Curve);
    table.push(Pi);
  }
  return table;
}

function random(Curve: AffineTwistedCurve) {
  let x = Curve.Field.random();
  return simpleMapToCurve(x, Curve);
}

/**
 * Given an x coordinate (base field element), increment it until we find one with
 * a y coordinate that satisfies the curve equation, and return the point.
 *
 * If the curve has a cofactor, multiply by it to get a point in the correct subgroup.
 */
function simpleMapToCurve(x: bigint, Curve: AffineTwistedCurve) {
  const F = Curve.Field;
  let y: bigint | undefined = undefined;

  // increment x until we find a y coordinate
  while (y === undefined) {
    x = F.add(x, 1n);
    // solve y^2 = (1 - a * x^2)/(1 - d * x^2)
    let x2 = F.square(x);
    let num = F.sub(1n, F.mul(x2, Curve.a));
    let den = F.sub(1n, F.mul(x2, Curve.d));
    if (den == 0n) continue;
    let y2 = F.div(num, den)!; // guaranteed that den has an inverse
    y = F.sqrt(y2);
  }

  let p = { x, y };

  // clear cofactor
  if (Curve.hasCofactor) {
    p = Curve.scale(p, Curve.cofactor!);
  }
  return p;
}

namespace Eddsa {
  /**
   * EdDSA signature consisting of a compressed curve point R and the scalar s.
   */
  export type Signature = { R: Field3; s: Field3 };
  export type signature = { R: bigint; s: bigint };
}

const EddsaSignature = {
  from({ R, s }: Eddsa.signature): Eddsa.Signature {
    return { R: Field3.from(R), s: Field3.from(s) };
  },
  toBigint({ R, s }: Eddsa.Signature): Eddsa.signature {
    return { R: Field3.toBigint(R), s: Field3.toBigint(s) };
  },
  isConstant: (S: Eddsa.Signature) => Provable.isConstant(EddsaSignature, S),

  /**
   * Parse an EdDSA signature from a raw 130-character hex string (64 bytes + "0x").
   */
  fromHex(rawSignature: string): Eddsa.Signature {
    // Validate input format
    let prefix = rawSignature.slice(0, 2);
    let signature = rawSignature.slice(2);
    if (prefix !== '0x' || signature.length !== 128) {
      throw new Error(
        `Signature.fromHex(): Invalid signature, expected hex string 0x... of length 130.`
      );
    }

    // Split the signature into R and s components
    const Rhex = signature.slice(0, 64); // First 32 bytes (64 hex chars for R)
    const Shex = signature.slice(64); // Last 32 bytes (64 hex chars for s)
    const R = BigInt(`0x${Rhex}`); // R value as a bigint
    const s = BigInt(`0x${Shex}`); // s value as a bigint

    if (s < 0 || s >= Curve.order) {
      throw new Error(`Invalid s value: must be a scalar modulo curve order.`);
    }

    Point.fromHex(Rhex); // Check that R represents a valid point

    return Eddsa.Signature.from({ R, s });
  },

  provable: provable({ R: Field3, s: Field3 }),
};

// Parameters used in Ed25519 (EdDSA algorithm for edwards25519 curve)
// https://datatracker.ietf.org/doc/html/rfc8032#section-5.1
const edwards25519Params: TwistedCurveParams = {
  name: 'edwards25519',
  modulus: (1n << 255n) - 19n, // 2^255 - 19
  order: 0x1000000000000000000000000000000014def9dea2f79cd65812631a5cf5d3edn, //2^252 + 27742317777372353535851937790883648493,
  cofactor: 8n,
  generator: {
    x: 0x216936d3cd6e53fec0a4e231fdd6dc5c692cc7609525a7b2c9562d608f25d51an, // <=> 15112221349535400772501151409588531511454012693041857206046113283949847762202
    y: 0x6666666666666666666666666666666666666666666666666666666666666658n, // <=> 4/5 mod p <=> 46316835694926478169428394003475163141307993866256225615783033603165251855960
  },
  a: 0x7fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffecn, // <=> -1 mod p <=> 57896044618658097711785492504343953926634992332820282019728792003956564819948
  d: 0x52036cee2b6ffe738cc740797779e89800700a4d4141d8ab75eb4dca135978a3n, // -121665/121666 mod p <=> 37095705934669439343138083508754565189542113879843219016388785533085940283555
};

const TwistedCurves = {
  Edwards25519: edwards25519Params,
};

/** EdDSA over Edwards25519 */
const Curve = createAffineTwistedCurve(edwards25519Params);
const basePoint = Point.from(Curve.one);

/**
 * Encode a point of the Edwards25519 curve into its compressed representation.
 *
 * @param input Point with {@link Field3} coordinates {x, y}
 * @returns 32-byte compressed representation of the point as {@link Field3}
 */
function encode(input: Point): Field3 {
  let p = Curve.Field.modulus;
  // https://www.rfc-editor.org/rfc/pdfrfc/rfc8032.txt.pdf Section 5.1.2
  let witnesses = Provable.witnessFields(8, () => {
    let x = Field3.toBigint(input.x);
    let y = Field3.toBigint(input.y);
    let x_lsb = x & 1n; // parity bit for x
    let x_masked = (x >> 1n) * 2n; // x with parity bit removed
    let y_msb = (y >> 255n) & 1n; // most significant bit of y
    let y_masked = y & ((1n << 255n) - 1n); // mask most significant bit

    return [x_lsb, ...Field3.split(x_masked), y_msb, ...Field3.split(y_masked)];
  });

  let [
    x_lsb,
    x_masked0,
    x_masked1,
    x_masked2,
    y_msb,
    y_masked0,
    y_masked1,
    y_masked2,
  ] = witnesses;

  x_lsb.assertBool('Parity bit of x coordinate is not a bit');
  y_msb.assertBool('MSB of y coordinate is not a bit');

  let x_masked: Field3 = [x_masked0, x_masked1, x_masked2];
  let y_masked: Field3 = [y_masked0, y_masked1, y_masked2];

  let x_lsb3: Field3 = [x_lsb, new Field(0n), new Field(0n)];
  let y_msb3: Field3 = [y_msb, new Field(0n), new Field(0n)];

  ForeignField.assertEquals(input.x, ForeignField.add(x_lsb3, x_masked, p));
  ForeignField.assertEquals(
    input.y,
    ForeignField.add(
      ForeignField.mul(y_msb3, Field3.from(1n << 255n), p),
      y_masked,
      p
    )
  );

  let enc = ForeignField.add(
    ForeignField.mul(x_lsb3, Field3.from(1n << 255n), p),
    y_masked,
    p
  );
  return enc;
}

/**
 * Decode a little-endian 32-byte compressed representation of Edwards25519 as
 * the point.
 *
 * @param compressed: 32-byte compressed representation of the point as {@link Field3}
 * @returns Point with {@link Field3} coordinates {x, y}
 */
function decode(input: UInt8[]): Point {
  if (input.length !== 32) {
    throw new Error(
      `Invalid compressed point: expected 32 bytes, got ${input.length}.`
    );
  }

  let p = Curve.modulus;

  // https://www.rfc-editor.org/rfc/pdfrfc/rfc8032.txt.pdf Section 5.1.3
  let witnesses = Provable.witnessFields(11, () => {
    let bytes = input.map((byte) => byte.toBigInt());
    // most significant byte of input is the parity bit for x
    const x_par = bytes[31] >> 7n;
    bytes[31] &= 0b01111111n;
    const y_msb = bytes[31];

    const y = bytes.reduce((acc, byte, i) => acc + (byte << BigInt(i * 8)), 0n);
    assert(y < p, 'Decoding failed: y coordinate larger than the field size');
    const x = recoverX(y, x_par);

    // to check parity bit of x
    const aux = (x - x_par) / 2n;

    return [
      ...Field3.split(aux),
      x_par,
      ...Field3.split(x),
      y_msb,
      ...Field3.split(y),
    ];
  });

  let [aux0, aux1, aux2, x_par, x0, x1, x2, y_msb, y0, y1, y2] = witnesses;
  let x_0limbs: Field3 = [x_par, Field.from(0n), Field.from(0n)];
  let aux: Field3 = [aux0, aux1, aux2];
  let x: Field3 = [x0, x1, x2];
  let y: Field3 = [y0, y1, y2];

  ForeignField.assertLessThan(y, p);
  // check x_0 is a bit
  x_par.assertBool('Parity bit of x coordinate is not a bit');

  // check y_msb shape
  input[31].value.assertEquals(y_msb.add(x_par.mul(128n)));

  // check y decomposition
  let input_ = input.slice();
  input_[31] = UInt8.from(y_msb);
  ForeignField.assertEquals(y, toField3(input_, p));

  // check (x, y) is on the curve
  TwistedCurve.assertOnCurve({ x, y }, Curve);

  // check parity/sign of x
  ForeignField.assertEquals(
    ForeignField.add(aux, aux, p),
    ForeignField.sub(x, x_0limbs, p)
  );

  // if x is zero, x_0 must be 0
  Provable.if(
    ForeignField.equals(x, 0n, p),
    Bool,
    x_par.equals(1n).not(),
    new Bool(true)
  );
  // check sign of x
  // if x_0 is 1, x must be odd (negative)
  Provable.if(x_par.equals(1n), Bool, x[0].equals(1n), new Bool(true));

  return { x, y };
}

/**
 * Generate a new EdDSA public key from a private key that is a random 32-byte
 * random seed.
 *
 * https://www.rfc-editor.org/rfc/pdfrfc/rfc8032.txt.pdf Section 5.1.5
 *
 * @param privateKey: 32-byte random seed
 * @returns the public key as a 32-byte encoded curve point,
 *          and the full SHA2-512 digest of the private key
 */
function keygenEddsa(privateKey: UInt8[]): [Field3, Bytes] {
  // TODO: use arrays instead of bigints?
  if (privateKey.length > 32) {
    throw new Error(`Invalid length of EdDSA private key: ${privateKey}.`);
  }
  // hash private key with SHA2-512
  const h = SHA2.hash(512, privateKey);
  // only need lowest 32 bytes to generate the public key
  let buffer = h.bytes.slice(0, 32);
  // prune buffer
  buffer[0] = UInt8.from(
    Gadgets.and(buffer[0].value, Field.from(0b11111000), 8)
  ); // clear lowest 3 bits
  buffer[31] = UInt8.from(
    Gadgets.and(buffer[31].value, Field.from(0b01111111), 8)
  ); // clear highest bit
  buffer[31] = UInt8.from(
    Gadgets.or(buffer[31].value, Field.from(0b01000000), 8)
  ); // set second highest bit

  // NOTE: despite clearing the top bit,
  //       the scalar could be larger than a native field element

  // read scalar from buffer (initially laid out as little endian)
  const f = Curve.Field.modulus;
  const s = toField3(buffer, f);

  return [encode(TwistedCurve.scale(s, basePoint, Curve)), h];
}

/**
 * Sign a message using Ed25519 (EdDSA over Edwards25519 curve).
 *
 * https://www.rfc-editor.org/rfc/pdfrfc/rfc8032.txt.pdf Section 5.1.6
 *
 * @param privateKey: 32-byte random seed
 * @param message: arbitrary length message to be signed
 * @returns the 64-bit signature composed by 32-bytes corresponding to a
 *          compressed curve point and a 32-byte scalar
 */
function signEddsa(
  privateKey: bigint,
  message: (bigint | number)[] | Uint8Array
): Eddsa.Signature {
  const L = Curve.order;
  let key = fromBigint(privateKey);
  const [publicKey, h] = keygenEddsa(key);
  // secret scalar obtained from first half of the digest
  const scalar = h.bytes.slice(0, 32);
  // prefix obtained from second half of the digest
  const prefix = h.bytes.slice(32, 64);

  // Hash the prefix concatenated with the message to obtain 64 bytes, that
  // need to be interpreted as little endian and reduced modulo the curve order
  const r = toField3(SHA2.hash(512, [...prefix, ...message]).bytes, L);

  // R is the encoding of the point resulting from [r]B
  let R = encode(TwistedCurve.scale(r, basePoint, Curve));

  // Hash the encoding concatenated with the public key and the message to
  // obtain a 64-byte digest that needs to be interpreted as little endian
  // and reduced modulo the curve order
  const k = toField3(
    SHA2.hash(512, [
      ...fromField3(R).flat(),
      ...fromField3(publicKey).flat(),
      ...message,
    ]).bytes,
    L
  );

  let s = ForeignField.add(r, ForeignField.mul(k, toField3(scalar, L), L), L);

  return { R, s };
}

function verifyEddsa(
  signature: Eddsa.Signature,
  message: UInt8[],
  publicKey: Field3
): Bool {
  let { R, s } = signature;

  let { x, y } = decode(fromField3(R));
  let A = decode(fromField3(publicKey));

  ForeignField.assertLessThanOrEqual(s, Curve.order);

  let k = SHA2.hash(512, [
    ...fromField3(R).flat(),
    ...fromField3(publicKey).flat(),
    ...message.flat(),
  ]).bytes;

  // Check [s]B = R + [k]A
  return Provable.equal(
    Point,
    TwistedCurve.scale(s, basePoint, Curve),
    TwistedCurve.add(
      { x, y },
      TwistedCurve.scale(toField3(k, Curve.Field.modulus), A, Curve),
      Curve
    )
  );
}

const Eddsa = {
  sign: signEddsa,
  verify: verifyEddsa,
  Signature: EddsaSignature,
};

// https://www.rfc-editor.org/rfc/pdfrfc/rfc8032.txt.pdf Section 5.1.3
function recoverX(y: bigint, x_0: bigint): bigint {
  const p = Curve.modulus;
  const u = y * y - 1n;
  const v = Curve.d * y * y - Curve.a;
  const candidate_x = (u * v) ^ (3n * ((u * v) ^ 7n)) ^ ((p - 5n) / 8n);

  let aux = mod((v * candidate_x) ^ 2n, p);

  let x =
    aux === u
      ? candidate_x
      : aux === -u
      ? (candidate_x * 2n) ^ ((p - 1n) / 4n)
      : (() => {
          throw new Error(
            `Decoding failed: no square root x exists for y value: ${y}.`
          );
        })();

  // Use the parity bit to select the correct sign for x
  if (x === 0n && x_0 === 1n) {
    throw new Error(`Invalid x value: x is zero but parity bit is 1.`);
  } else if (x % 2n !== x_0) {
    x = p - x;
  }

  return x;
}

class ForeignTwisted {
  x: AlmostForeignField;
  y: AlmostForeignField;

  /**
   * Create a new {@link ForeignTwisted} from an object representing the (affine
   * twisted) x and y coordinates.
   *
   * Note: Inputs must be range checked if they originate from a different field
   * with a different modulus or if they are not constants. Please refer to the
   * {@link ForeignField} constructor comments for more details.
   *
   * @example
   * ```ts
   * let x = new ForeignTwisted({ x: 0n, y: 1n });
   * ```
   *
   * **Warning**: This fails for a constant input which does not represent an
   *              actual point on the curve or not in the subgroup.
   *
   * **Note**: For now, only the edwards25519 curve is supported.
   */
  constructor(g: {
    x: AlmostForeignField | Field3 | bigint | number;
    y: AlmostForeignField | Field3 | bigint | number;
  }) {
    this.x = new this.Constructor.Field(g.x);
    this.y = new this.Constructor.Field(g.y);
    // don't allow constants that aren't on the curve or in the prime subgroup
    if (this.isConstant()) {
      this.assertOnCurve();
      this.assertInSubgroup();
    }
  }

  /**
   * Coerce the input to a {@link ForeignTwisted}.
   */
  static from(g: ForeignTwisted | FlexiblePoint) {
    if (g instanceof this) return g;
    return new this(g);
  }

  /**
   * Parses a hexadecimal string representing an uncompressed elliptic curve point
   * and coerces it into a {@link ForeignTwisted} point using big-endian byte order.
   *
   * The method extracts the x and y coordinates from the provided hex string and
   * verifies that the resulting point lies on the curve.
   *
   * **Note:** This method only supports uncompressed elliptic curve points, which
   * are 65 bytes in total (1-byte prefix + 32 bytes for x + 32 bytes for y).
   *
   * @param hex - The hexadecimal string representing the uncompressed elliptic curve point.
   * @returns - A point on the foreign curve, parsed from the given hexadecimal string.
   *
   * @throws - Throws an error if the input is not a valid public key.
   *
   * @example
   * ```ts
   * class Edwards25519 extends createForeignCurveTwisted(Crypto.TwistedCurveParams.Edwards25519) {}
   *
   * // Example hex string for uncompressed point
   * const publicKeyHex = '04f8b8db25c619d0c66b2dc9e97ecbafafae...';
   * const point = Edwards25519.fromHex(publicKeyHex);
   * ```
   *
   * **Important:** This method is only designed to handle uncompressed elliptic curve points in hex format.
   */
  static fromHex(hex: string) {
    // trim the '0x' prefix if present
    if (hex.startsWith('0x')) {
      hex = hex.slice(2);
    }

    const bytes = Bytes.fromHex(hex).toBytes();
    const sizeInBytes = Math.ceil(this.Bigint.Field.sizeInBits / 8);

    // extract x and y coordinates from the byte array
    const xBytes = bytes.subarray(0, sizeInBytes); // first `sizeInBytes` bytes for x-coordinate
    const yBytes = bytes.subarray(sizeInBytes, 2 * sizeInBytes); // next `sizeInBytes` bytes for y-coordinate

    // convert byte arrays to bigint
    const x = BigInt('0x' + Bytes.from(xBytes).toHex());
    const y = BigInt('0x' + Bytes.from(yBytes).toHex());

    // construct the point on the curve using the x and y coordinates
    let P = this.from({ x, y });

    // ensure that the point is on the curve
    P.assertOnCurve();

    return P;
  }

  /**
   * The constant generator point.
   */
  static get generator() {
    return new this(this.Bigint.one);
  }
  /**
   * The size of the curve's base field.
   */
  static get modulus() {
    return this.Bigint.modulus;
  }
  /**
   * The size of the curve's base field.
   */
  get modulus() {
    return this.Constructor.Bigint.modulus;
  }

  /**
   * Checks whether this curve point is constant.
   *
   * See {@link FieldVar} to understand constants vs variables.
   */
  isConstant() {
    return Provable.isConstant(this.Constructor, this);
  }

  /**
   * Convert this curve point to a point with bigint coordinates.
   */
  toBigint() {
    return this.Constructor.Bigint.from({
      x: this.x.toBigInt(),
      y: this.y.toBigInt(),
    });
  }

  /**
   * Twisted elliptic curve addition (complete)
   *
   * ```ts
   * let r = p.add(q); // r = p + q
   * ```
   */
  add(h: ForeignTwisted | FlexiblePoint) {
    let Curve = this.Constructor.Bigint;
    let h_ = this.Constructor.from(h);
    let p = TwistedCurve.add(toPoint(this), toPoint(h_), Curve);
    return new this.Constructor(p);
  }

  /**
   * Twisted elliptic curve doubling.
   *
   * @example
   * ```ts
   * let r = p.double(); // r = 2 * p
   * ```
   */
  double() {
    let Curve = this.Constructor.Bigint;
    let p = TwistedCurve.double(toPoint(this), Curve);
    return new this.Constructor(p);
  }

  /**
   * Twisted elliptic curve negation.
   *
   * @example
   * ```ts
   * let r = p.negate(); // r = -p
   * ```
   */
  negate(): ForeignTwisted {
    return new this.Constructor({ x: this.x.neg(), y: this.y });
  }

  /**
   * Twisted elliptic curve scalar multiplication, where the scalar is represented as a {@link ForeignField} element.
   *
   * @example
   * ```ts
   * let r = p.scale(s); // r = s * p
   * ```
   */
  scale(scalar: AlmostForeignField | bigint | number) {
    let Curve = this.Constructor.Bigint;
    let scalar_ = this.Constructor.Scalar.from(scalar);
    let p = TwistedCurve.scale(scalar_.value, toPoint(this), Curve);
    return new this.Constructor(p);
  }

  static assertOnCurve(g: ForeignTwisted) {
    TwistedCurve.assertOnCurve(toPoint(g), this.Bigint);
  }

  /**
   * Assert that this point lies on the elliptic curve, which means it satisfies the equation
   * `a * x^2 + y^2 = 1 + d * x^2 * y^2`
   */
  assertOnCurve() {
    this.Constructor.assertOnCurve(this);
  }

  static assertInSubgroup(g: ForeignTwisted) {
    if (this.Bigint.hasCofactor) {
      TwistedCurve.assertInSubgroup(toPoint(g), this.Bigint);
    }
  }

  /**
   * Assert that this point lies in the prime subgroup, which means that scaling
   * the point by the curve order results in a nonzero point.
   */
  assertInSubgroup() {
    this.Constructor.assertInSubgroup(this);
  }

  /**
   * Check that this is a valid element of the target subgroup of the curve:
   * - Check that the coordinates are valid field elements
   * - Use {@link assertOnCurve()} to check that the point lies on the curve
   * - If the curve has cofactor unequal to 1, use {@link assertInSubgroup()}.
   */
  static check(g: ForeignTwistedNotNeeded) {
    multiRangeCheck(g.x.value);
    multiRangeCheck(g.y.value);
    this.assertOnCurve(g);
    this.assertInSubgroup(g);
  }

  // dynamic subclassing infra
  get Constructor() {
    return this.constructor as typeof ForeignTwisted;
  }
  static _Bigint?: AffineTwistedCurve;
  static _Field?: typeof AlmostForeignField;
  static _Scalar?: typeof AlmostForeignField;
  static _provable?: ProvablePureExtended<
    ForeignTwisted,
    { x: bigint; y: bigint },
    { x: string; y: string }
  >;

  /**
   * Curve arithmetic on JS bigints.
   */
  static get Bigint() {
    assert(this._Bigint !== undefined, 'ForeignTwisted not initialized');
    return this._Bigint;
  }
  /**
   * The base field of this curve as a {@link ForeignField}.
   */
  static get Field() {
    assert(this._Field !== undefined, 'ForeignTwisted not initialized');
    return this._Field;
  }
  /**
   * The scalar field of this curve as a {@link ForeignField}.
   */
  static get Scalar() {
    assert(this._Scalar !== undefined, 'ForeignTwisted not initialized');
    return this._Scalar;
  }
  /**
   * `Provable<ForeignCurve>`
   */
  static get provable() {
    assert(this._provable !== undefined, 'ForeignTwisted not initialized');
    return this._provable;
  }
}

class ForeignTwistedNotNeeded extends ForeignTwisted {
  constructor(g: {
    x: AlmostForeignField | Field3 | bigint | number;
    y: AlmostForeignField | Field3 | bigint | number;
  }) {
    super(g);
  }

  static check(g: ForeignTwistedNotNeeded) {
    multiRangeCheck(g.x.value);
    multiRangeCheck(g.y.value);
    this.assertOnCurve(g);
    this.assertInSubgroup(g);
  }
}

/**
 * Create a class representing a twisted elliptic curve group, which is different
 * from the native {@link Group}.
 *
 * ```ts
 * const Curve = createForeignTwisted(Crypto.TwistedCurveParams.Edwards25519);
 * ```
 *
 * `createForeignTwisted(params)` takes curve parameters {@link TwistedCurveParams} as input.
 * We support `modulus` and `order` to be prime numbers up to 259 bits.
 *
 * The returned {@link ForeignTwistedNotNeeded} class represents a curve point
 * (including zero) and supports standard elliptic curve operations like point
 * addition and scalar multiplication.
 *
 * {@link ForeignTwistedNotNeeded} also includes to associated foreign fields:
 * `ForeignCurve.Field` and `ForeignCurve.Scalar`, see {@link createForeignField}.
 */
function createForeignTwisted(
  params: TwistedCurveParams
): typeof ForeignTwisted {
  assert(
    params.modulus > l2Mask + 1n,
    'Base field moduli smaller than 2^176 are not supported'
  );

  const FieldUnreduced = createForeignField(params.modulus);
  const ScalarUnreduced = createForeignField(params.order);
  class Field extends FieldUnreduced.AlmostReduced {}
  class Scalar extends ScalarUnreduced.AlmostReduced {}

  const BigintCurve = createAffineTwistedCurve(params);

  class Curve extends ForeignTwisted {
    static _Bigint = BigintCurve;
    static _Field = Field;
    static _Scalar = Scalar;
    static _provable = provableFromClass(Curve, { x: Field, y: Field });
  }

  return Curve;
}

function toField3(bytes: UInt8[], mod: bigint): Field3 {
  // TODO: more efficient implementation
  assert(mod < 1n << 259n, 'Foreign modulus must fit in 259 bits');
  return bytes
    .slice() // copy the array to prevent mutation
    .reverse()
    .map((b) => [Field.from(b.value), Field.from(0n), Field.from(0n)] as Field3)
    .reduce((acc, byte) =>
      ForeignField.add(
        ForeignField.mul(Field3.from(acc), Field3.from(256n), mod),
        Field3.from(byte),
        mod
      )
    );
}

function fromBigint(x: bigint, bytelength: number = 32): UInt8[] {
  assert(x < 1n << BigInt(bytelength * 8), 'Input does not fit in bytelength');
  let bytes = Array.from(
    { length: bytelength },
    (_, k) => new UInt8((x >> BigInt(8 * k)) & 0xffn)
  );
  return bytes;
}

function fromField3(x: Field3): UInt8[] {
  const limbBytes = Number(Gadgets.Constants.l) / 8;
  return [
    fromField(x[0], limbBytes),
    fromField(x[1], limbBytes),
    fromField(x[2], limbBytes - 1), // only 256 bits
  ].flat();
}
/**
 * Returns an array of {@link UInt8} elements representing this field element
 * as little endian ordered bytes.
 *
 * If the optional `bytelength` argument is used, it proves that the field
 * element fits in `bytelength` bytes. The length has to be between 0 and 32,
 * and the method throws if it isn't.
 *
 * **Warning**: The cost of this operation in a zk proof depends on the
 * `bytelength` you specify, which by default is 32 bytes. Prefer to pass a
 * smaller `bytelength` if possible.
 *
 * @param input - the field element to convert to bytes.
 * @param bytelength - the number of bytes to fit the element. If the element
 *                     does not fit in `length` bits, the functions throws an
 *                     error.
 *
 * @return An array of {@link UInt8} element representing this {@link Field} in
 *         little endian encoding.
 */
function fromField(input: Field, bytelength: number = 32): UInt8[] {
  checkBitLength('Field.toBytes()', bytelength, 32 * 8);
  if (input.isConstant()) {
    if (input.toBigInt() >= 1n << (BigInt(bytelength) * 8n)) {
      throw Error(`toOctets(): ${input} does not fit in ${bytelength} bytes`);
    }
    let x = input.toBigInt();
    return Array.from(
      { length: bytelength },
      (_, k) => new UInt8((x >> BigInt(8 * k)) & 0xffn)
    );
  }
  let bytes = Provable.witness(Provable.Array(UInt8, bytelength), () => {
    let x = input.toBigInt();
    return Array.from(
      { length: bytelength },
      (_, k) => new UInt8((x >> BigInt(8 * k)) & 0xffn)
    );
  });
  let field = bytes
    .reverse()
    .map((x) => x.value)
    .reduce((acc, byte) => acc.mul(256).add(byte));

  field.assertEquals(
    input,
    `toOctets(): incorrect decomposition into ${bytelength} bytes`
  );
  return bytes;
}

function checkBitLength(name: string, length: number, maxLength: number) {
  if (length > maxLength)
    throw Error(
      `${name}: bit length must be ${maxLength} or less, got ${length}`
    );
  if (length < 0)
    throw Error(`${name}: bit length must be non-negative, got ${length}`);
}

function assertPositiveInteger(n: number, message: string) {
  if (!Number.isInteger(n) || n <= 0) throw Error(message);
}
