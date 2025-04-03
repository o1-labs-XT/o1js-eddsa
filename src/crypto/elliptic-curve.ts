import { assert, Core } from 'o1js';

export {
  AffineTwistedCurve,
  GroupAffineTwisted,
  affineTwistedAdd,
  affineTwistedDouble,
  affineTwistedZero,
  createAffineTwistedCurve,
  TwistedCurveParams,
};

const { mod, inverse, createField, p } = Core.FiniteField;
const { bigIntToBits } = Core.BigIntHelpers;
type GroupProjective = Core.EllipticCurve.GroupProjective;

type GroupAffineTwisted = { x: bigint; y: bigint };

const projectiveTwistedZero = { x: 0n, y: 1n, z: 1n };
const affineTwistedZero: GroupAffineTwisted = {
  x: 0n,
  y: 1n,
};

/**
 * Parameters defining an elliptic curve in twisted Edwards form with affine coordinates
 * ax^2 + y^2 = 1 + dx^2y^2
 */
type TwistedCurveParams = {
  /**
   * Human-friendly name for the curve
   */
  name: string;
  /**
   * Base field modulus
   */
  modulus: bigint;
  /**
   * Scalar field modulus = group order
   */
  order: bigint;
  /**
   * Cofactor = size of EC / order
   *
   * This can be left undefined if the cofactor is 1.
   */
  cofactor?: bigint;
  /**
   * Generator point
   */
  generator: { x: bigint; y: bigint };
  /**
   * The `a` parameter in the curve equation ax^2 + y^2 = 1 + dx^2y^2
   */
  a: bigint;
  /**
   * The `d` parameter in the curve equation ax^2 + y^2 = 1 + dx^2y^2
   */
  d: bigint;

  endoBase?: bigint;
  endoScalar?: bigint;
};

function affineTwistedOnCurve(
  g: GroupAffineTwisted,
  p: bigint,
  a: bigint,
  d: bigint
) {
  if (affineTwistedIsZero(g, p)) return true;
  // a * x^2 + y^2 = 1 + d * x^2 * y^2
  const { x, y } = g;
  let x2 = mod(x * x, p);
  let y2 = mod(y * y, p);
  return mod(a * x2 + y2 - 1n - d * x2 * y2, p) === 0n;
}

// https://www.hyperelliptic.org/EFD/g1p/auto-twisted.html
function affineTwistedAdd(
  g: GroupAffineTwisted,
  h: GroupAffineTwisted,
  p: bigint,
  a: bigint,
  d: bigint
): GroupAffineTwisted {
  if (affineTwistedIsZero(g, p)) return h;
  if (affineTwistedIsZero(h, p)) return g;

  let { x: x1, y: y1 } = g;
  let { x: x2, y: y2 } = h;

  if (y1 === y2) {
    // g + g --> double
    if (x1 === x2) return affineTwistedDouble(g, p, a, d);
    // g - g --> return zero
    if (x1 === mod(p - x2, p)) return affineTwistedZero;
  }

  // x3 = (x1 * y2 + y1 * x2) / (1 + d * x1 * x2 * y1 * y2)
  // y3 = (y1 * y2 - a * x1 * x2) / (1 - d * x1 * x2 * y1 * y2)
  let x1x2 = mod(x1 * x2, p);
  let y1y2 = mod(y1 * y2, p);
  let x1y2 = mod(x1 * y2, p);
  let y1x2 = mod(y1 * x2, p);
  let ax1x2 = mod(a * x1x2, p);

  let x3Num = mod(x1y2 + y1x2, p);
  let y3Num = mod(y1y2 - ax1x2, p);

  let dx1x2y1y2 = mod(d * x1x2 * y1y2, p);

  let x3Denom = inverse(mod(1n + dx1x2y1y2, p), p);
  if (x3Denom === undefined)
    throw Error('X denominator used in twisted addition is 0');

  let y3Denom = inverse(mod(1n - dx1x2y1y2, p), p);
  if (y3Denom === undefined)
    throw Error('Y denominator used in twisted addition is 0');

  let x3 = mod(x3Num * x3Denom, p);
  let y3 = mod(y3Num * y3Denom, p);

  return { x: x3, y: y3 };
}

// https://www.hyperelliptic.org/EFD/g1p/auto-twisted.html
function affineTwistedDouble(
  g: GroupAffineTwisted,
  p: bigint,
  a: bigint,
  d: bigint
): GroupAffineTwisted {
  let { x: x1, y: y1 } = g;

  if (g == affineTwistedZero) return g;

  // x3 = 2*x1*y1 / (1 + d * x1^2 * y1^2)
  // y3 = (y1^2 - a * x1^2) / (1 - d * x1^2 * y1^2)
  let x1x1 = x1 * x1;
  let y1y1 = y1 * y1;
  let x1y1 = x1 * y1;

  let x3Num = mod(2n * x1y1, p);
  let y3Num = mod(y1y1 - a * x1x1, p);

  let dx1x1y1y1 = mod(d * x1x1 * y1y1, p);

  let x3Den = inverse(1n + dx1x1y1y1, p);
  if (x3Den === undefined) throw Error('impossible');
  let y3Den = inverse(1n - dx1x1y1y1, p);
  if (y3Den === undefined) throw Error('impossible');

  let x3 = mod(x3Num * x3Den, p);
  let y3 = mod(y3Num * y3Den, p);

  return { x: x3, y: y3 };
}

function affineTwistedNegate(
  g: GroupAffineTwisted,
  p: bigint
): GroupAffineTwisted {
  if (g == affineTwistedZero) return g;
  return { x: g.x === 0n ? 0n : p - g.x, y: g.y };
}

function affineTwistedScale(
  g: GroupAffineTwisted,
  s: bigint | boolean[],
  p: bigint,
  a: bigint,
  d: bigint
) {
  let gProj = projectiveFromAffineTwisted(g);
  let sgProj = projectiveTwistedScale(gProj, s, p, a, d);
  return projectiveToAffineTwisted(sgProj, p);
}

// https://www.hyperelliptic.org/EFD/g1p/auto-twisted-projective.html
// https://eprint.iacr.org/2008/013.pdf Section 6
function projectiveTwistedAdd(
  g: GroupProjective,
  h: GroupProjective,
  p: bigint,
  a: bigint,
  d: bigint
): GroupProjective {
  let { x: X1, y: Y1, z: Z1 } = g;
  let { x: X2, y: Y2, z: Z2 } = h;

  // A = Z1 * Z2
  let A = mod(Z1 * Z2, p);
  // B = A^2
  let B = mod(A * A, p);
  // C = X1 * X2
  let C = mod(X1 * X2, p);
  // D = Y1 * Y2
  let D = mod(Y1 * Y2, p);
  // E = d * C * D
  let E = mod(d * C * D, p);
  // F = B - E
  let F = mod(B - E, p);
  // G = B + E
  let G = mod(B + E, p);

  return {
    x: mod(A * F * ((X1 + Y1) * (X2 + Y2) - C - D), p),
    y: mod(A * G * (D - a * C), p),
    z: mod(F * G, p),
  };
}

// https://www.hyperelliptic.org/EFD/g1p/auto-twisted-projective.html
// https://eprint.iacr.org/2008/013.pdf Section 6
function projectiveTwistedDouble(
  g: GroupProjective,
  p: bigint,
  a: bigint
): GroupProjective {
  let { x: X, y: Y, z: Z } = g;

  // B = (X + Y)^2
  let B = mod((X + Y) ** 2n, p);
  // C = X^2
  let C = mod(X * X, p);
  // D = Y^2
  let D = mod(Y * Y, p);
  // E = a * C
  let E = mod(a * C, p);
  // F = E + D
  let F = mod(E + D, p);
  // H = Z^2
  let H = mod(Z * Z, p);
  // J =  F - 2 * H
  let J = mod(F - 2n * H, p);

  return {
    x: mod((B - C - D) * J, p),
    y: mod(F * (E - D), p),
    z: mod(F * J, p),
  };
}

function projectiveTwistedScale(
  g: GroupProjective,
  x: bigint | boolean[],
  p: bigint,
  a: bigint,
  d: bigint
) {
  let bits = typeof x === 'bigint' ? bigIntToBits(x) : x;
  let h = projectiveTwistedZero;
  for (let bit of bits) {
    if (bit) h = projectiveTwistedAdd(h, g, p, a, d);
    g = projectiveTwistedDouble(g, p, a);
  }
  return h;
}

function projectiveFromAffineTwisted(g: GroupAffineTwisted): GroupProjective {
  if (affineTwistedIsZero(g, p)) return projectiveTwistedZero;
  return { x: g.x, y: g.y, z: 1n };
}

// The affine twisted curve with equation
// a * x^2 + y^2 = 1 + d * x^2 * y^2
// in projective coordinates is represented as
// a * X^2 * Z^2 + Y^2 Z^2 = Z^4 + d * X^2 * Y^2
// where x = X/Z, y = Y/Z, and Z ≠ 0
function projectiveToAffineTwisted(
  g: GroupProjective,
  p: bigint
): GroupAffineTwisted {
  let z = g.z;
  assert(
    z !== 0n,
    'projectiveToAffineTwisted: degenerate case, z must not be zero'
  );
  if (z === 1n) {
    // any other normalized affine form
    return { x: g.x, y: g.y };
  } else {
    let zinv = inverse(z, p)!; // we checked for z === 0, so inverse exists
    // x/z
    let x = mod(g.x * zinv, p);
    // y/z
    let y = mod(g.y * zinv, p);
    return { x, y };
  }
}

type AffineTwistedCurve = ReturnType<typeof createAffineTwistedCurve>;

function affineTwistedIsZero(g: GroupAffineTwisted, p: bigint): boolean {
  return mod(g.x, p) === 0n && mod(g.y, p) === 1n;
}

// auxiliary function to compute modular exponentiation (runs in O(log(exp)))
function modPow(base: bigint, exp: bigint, mod: bigint) {
  let result = 1n;
  base = base % mod; // Reduce base initially

  while (exp > 0n) {
    if (exp % 2n === 1n) {
      // If exponent is odd, multiply result
      result = (result * base) % mod;
    }
    base = (base * base) % mod; // Square the base
    exp = exp / 2n; // Halve the exponent
  }

  return result;
}

/** Creates twisted Edwards curves in affine cordinates of the form
 * a * x^2 + y^2 = 1 + d * x^2 * y^2
 * with a ≠ 0, d ≠ 0 and a ≠ d
 *
 * Warning: must be used only for curves without endomorphism, like edwards25519
 */
function createAffineTwistedCurve({
  name,
  modulus: p,
  order,
  cofactor,
  generator,
  a,
  d,
}: TwistedCurveParams) {
  let hasCofactor = cofactor !== undefined && cofactor !== 1n;

  const Field = createField(p);
  const Scalar = createField(order);
  const one = { ...generator };
  const Endo = undefined; // for edwards25519

  assert(a !== 0n, 'a must not be zero');
  assert(d !== 0n, 'd must not be zero');
  assert(a !== d, 'a must not be equal to d');
  // Euler's criterion: d is square iff d^((p - 1) / 2) = 1 mod p
  assert(modPow(d, (p - 1n) / 2n, p) != 1n, 'd must not be a square');

  return {
    name,
    /**
     * Arithmetic over the base field
     */
    Field,
    /**
     * Arithmetic over the scalar field
     */
    Scalar,

    modulus: p,
    order,
    a,
    d,
    cofactor,
    hasCofactor,

    zero: affineTwistedZero,
    one,

    hasEndomorphism: Endo !== undefined,
    get Endo() {
      if (Endo === undefined) throw Error(`no endomorphism defined on ${name}`);
      return Endo;
    },

    // Obtain a point from its affine representation, included the zero point.
    // NOTE: The method does not check that the point is on the curve nor the subgroup.
    from(g: { x: bigint; y: bigint }): GroupAffineTwisted {
      if (affineTwistedIsZero(g, p)) return affineTwistedZero;
      return { ...g };
    },

    equal(g: GroupAffineTwisted, h: GroupAffineTwisted) {
      if (affineTwistedIsZero(g, p) && affineTwistedIsZero(h, p)) {
        return true;
      } else if (affineTwistedIsZero(g, p) || affineTwistedIsZero(h, p)) {
        return false;
      } else {
        return mod(g.x - h.x, p) === 0n && mod(g.y - h.y, p) === 0n;
      }
    },
    isOnCurve(g: GroupAffineTwisted) {
      return affineTwistedOnCurve(g, p, a, d);
    },
    isInSubgroup(g: GroupAffineTwisted) {
      return affineTwistedIsZero(affineTwistedScale(g, order, p, a, d), p);
    },
    add(g: GroupAffineTwisted, h: GroupAffineTwisted) {
      return affineTwistedAdd(g, h, p, a, d);
    },
    double(g: GroupAffineTwisted) {
      return affineTwistedDouble(g, p, a, d);
    },
    negate(g: GroupAffineTwisted) {
      return affineTwistedNegate(g, p);
    },
    sub(g: GroupAffineTwisted, h: GroupAffineTwisted) {
      return affineTwistedAdd(g, affineTwistedNegate(h, p), p, a, d);
    },
    scale(g: GroupAffineTwisted, s: bigint | boolean[]) {
      return affineTwistedScale(g, s, p, a, d);
    },
    isZero(g: GroupAffineTwisted) {
      return affineTwistedIsZero(g, p);
    },
  };
}
