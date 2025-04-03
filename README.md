# eddsa-o1js

> ⚠️ **NOTE**: The current implementation produces circuits that are too large to fit within existing constraint limits. PRs that optimize this library are welcome!

![License: Apache-2.0](https://img.shields.io/badge/License-Apache%202.0-blue.svg)
[![npm version](https://img.shields.io/npm/v/eddsa-o1js.svg)](https://www.npmjs.com/package/eddsa-o1js)

A provable EdDSA signature verification library for [o1js](https://github.com/o1-labs/o1js), enabling zkApp developers to verify EdDSA signatures inside zk-SNARKs.

## Features

- **EdDSA Verification in ZK**: Verify EdDSA signatures within provable o1js code
- **Edwards25519 Support**: Built-in support for the popular Edwards25519 curve
- **Twisted Edwards Curves**: Implementation of twisted Edwards curves for efficient elliptic curve operations
- **Fully Composable**: Designed to work seamlessly with the o1js ecosystem

## Installation

```bash
npm install eddsa-o1js
```

## Usage

Here's a quick example of how to use eddsa-o1js to verify an EdDSA signature:

```typescript
import { ZkProgram, Bool, Bytes } from 'o1js';
import { createEddsa, createForeignTwisted, TwistedCurves } from 'eddsa-o1js';

// Create a custom Edwards25519 curve class
class Edwards25519 extends createForeignTwisted(TwistedCurves.Edwards25519) {}
class Scalar extends Edwards25519.Scalar {}
class Eddsa extends createEddsa(Edwards25519) {}
class Bytes32 extends Bytes(32) {}

// Define a ZkProgram that verifies EdDSA signatures
const eddsa = ZkProgram({
  name: 'eddsa',
  publicInput: Bytes32,
  publicOutput: Bool,

  methods: {
    verifyEddsa: {
      privateInputs: [Eddsa, Edwards25519],
      async method(
        message: Bytes32,
        signature: Eddsa,
        publicKey: Edwards25519
      ) {
        return {
          publicOutput: signature.verify(message, publicKey),
        };
      },
    },
  },
});

// Example: Generate a signature and verify it
async function run() {
  // Generate a keypair
  let privateKey = Edwards25519.Scalar.random();
  let publicKey = Edwards25519.generator.scale(privateKey);

  // Sign a message
  let message = Bytes32.fromString('Hello, o1js!');
  let signature = Eddsa.sign(message.toBytes(), privateKey.toBigInt());

  // Compile the program
  await eddsa.compile();

  // Verify the signature in zk
  let { proof } = await eddsa.verifyEddsa(message, signature, publicKey);

  // Check the result
  proof.publicOutput.assertTrue('signature verifies');
}
```

## Advanced Usage

For more detailed examples, please check the [examples directory](./examples):

- [Basic EdDSA signature verification](./examples/eddsa.ts)
- [Running a complete example](./examples/run.ts)

## API Reference

### Core Components

- `createEddsa(TwistedCurve)`: Factory function that creates an EdDSA implementation for a specific curve
- `createForeignTwisted(TwistedCurveParams)`: Creates a provable twisted Edwards curve implementation
- `TwistedCurves`: Contains parameters for common twisted Edwards curves (e.g., Edwards25519)

### Usage Patterns

1. Define your curve by extending the base implementation
2. Create your EdDSA implementation for that curve
3. Use the signature verification methods in your ZkProgram

## Development

```bash
# Build the project
npm run build
```

## Related Projects

- [o1js](https://github.com/o1-labs/o1js) - The main framework for writing zero-knowledge applications

## License

This project is licensed under the Apache License 2.0 - see the [LICENSE](./LICENSE) file for details.
