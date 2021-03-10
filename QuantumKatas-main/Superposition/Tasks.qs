// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT license.

namespace Quantum.Kata.Superposition {

    open Microsoft.Quantum.Diagnostics;
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Convert;
    open Microsoft.Quantum.Math;


    //////////////////////////////////////////////////////////////////
    // Welcome!
    //////////////////////////////////////////////////////////////////

    // "Superposition" quantum kata is a series of exercises designed
    // to get you familiar with programming in Q#.
    // It covers the following topics:
    //  - basic single-qubit and multi-qubit gates,
    //  - superposition,
    //  - flow control and recursion in Q#.

    // Each task is wrapped in one operation preceded by the description of the task.
    // Each task (except tasks in which you have to write a test) has a unit test associated with it,
    // which initially fails. Your goal is to fill in the blank (marked with // ... comment)
    // with some Q# code to make the failing test pass.

    // The tasks are given in approximate order of increasing difficulty; harder ones are marked with asterisks.

    //////////////////////////////////////////////////////////////////
    // Part I. Simple Gates
    //////////////////////////////////////////////////////////////////

    // Task 1.1. Plus state
    // Input: a qubit in the |0⟩ state.
    // Goal: prepare a |+⟩ state on this qubit (|+⟩ = (|0⟩ + |1⟩) / sqrt(2)).
    operation PlusState (q : Qubit) : Unit {
        // Hadamard gate H will convert |0⟩ state to |+⟩ state.
        // Type the following: H(q);
        // Then rebuild the project and rerun the tests - T01_PlusState should now pass!
        H(q);
    }


    // Task 1.2. Minus state
    // Input: a qubit in the |0⟩ state.
    // Goal: prepare a |-⟩ state on this qubit (|-⟩ = (|0⟩ - |1⟩) / sqrt(2)).
    operation MinusState (q : Qubit) : Unit {
        // In this task, as well as in all subsequent ones, you have to come up with the solution yourself.
        X(q);
        H(q);
    }


    // Task 1.3. Superposition of all basis vectors on two qubits
    // Input: two qubits in |00⟩ state (stored in an array of length 2).
    // Goal:  create the following state on these qubits: (|00⟩ + |01⟩ + |10⟩ + |11⟩) / 2.
    operation AllBasisVectors_TwoQubits (qs : Qubit[]) : Unit {
        // The following lines enforce the constraints on the input that you are given.
        // You don't need to modify them. Feel free to remove them, this won't cause your code to fail.
        Fact(Length(qs) == 2, "The array should have exactly 2 qubits.");
        H(qs[0]);
        H(qs[1]);
    }


    // Task 1.4. Superposition of basis vectors with phase flip.
    // Input: two qubits in |00⟩ state (stored in an array of length 2).
    // Goal:  create the following state on these qubits: (|00⟩ + |01⟩ + |10⟩ - |11⟩) / 2.
    operation AllBasisVectorWithPhaseFlip_TwoQubits (qs : Qubit[]) : Unit {
        //       H(qs[0]);
        //       H(qs[1]);
        AllBasisVectors_TwoQubits(qs);
        Controlled Z ([qs[0]], qs[1]);  // Just to apply it to |11>
    }


    // Task 1.5. Superposition of basis vectors with phases
    // Input: two qubits in |00⟩ state (stored in an array of length 2).
    // Goal:  create the following state on these qubits: (|00⟩ + i*|01⟩ - |10⟩ - i*|11⟩) / 2.
    operation AllBasisVectorsWithPhases_TwoQubits (qs : Qubit[]) : Unit {
        // 1/2 means H() is there
        H(qs[0]);
        H(qs[1]);
        // Need to apply -1 to |01> and -i to |10>
        Y(qs[0]);
        S(qs[1]);
    }


    // Task 1.6. Bell state
    // Input: two qubits in |00⟩ state (stored in an array of length 2).
    // Goal: create a Bell state |Φ⁺⟩ = (|00⟩ + |11⟩) / sqrt(2) on these qubits.
    operation BellState (qs : Qubit[]) : Unit {
        H(qs[0]);
        CNOT(qs[0], qs[1]);
    }


    // Task 1.7. All Bell states
    // Inputs:
    //      1) two qubits in |00⟩ state (stored in an array of length 2)
    //      2) an integer index
    // Goal: create one of the Bell states based on the value of index:
    //       0: |Φ⁺⟩ = (|00⟩ + |11⟩) / sqrt(2)
    //       1: |Φ⁻⟩ = (|00⟩ - |11⟩) / sqrt(2)
    //       2: |Ψ⁺⟩ = (|01⟩ + |10⟩) / sqrt(2)
    //       3: |Ψ⁻⟩ = (|01⟩ - |10⟩) / sqrt(2)
    operation AllBellStates (qs : Qubit[], index : Int) : Unit {
        // Case 0
        // Need alway to apply Hadamard
        H(qs[0]);
        
        // Start conditionals
        
        // Case 1
        if (index == 1) {
            Z(qs[0]);
        }
        // Case 2
        if (index == 2) {
            X(qs[1]);
        }
        
        // Case 3
        if (index == 3) {
            Z(qs[0]);
            X(qs[1]);
        }
        
        // Need to apply CNOT
        CNOT(qs[0], qs[1]);
    }


    // Task 1.8. Greenberger–Horne–Zeilinger state
    // Input: N qubits in |0...0⟩ state.
    // Goal: create a GHZ state (|0...0⟩ + |1...1⟩) / sqrt(2) on these qubits.
    operation GHZ_State (qs : Qubit[]) : Unit {
        // Hint: N can be found as Length(qs).
        //   |q1>------H----o------o---
        //                  |      |
        //   <q2|-----------+------|---
        //                         |
        //   |qn>------------------+---
        //
        
        // Apply H to q1
        H(qs[0]);

        // Loop through all other qubits
        //    for index in 1 .. (Length(qs)-1){
        //        CNOT(qs[0], qs[index]);
        //    }

        // Loop through all other qubits
        for q in Rest(qs) {
            CNOT(qs[0], q);
        }
    }


    // Task 1.9. Superposition of all basis vectors
    // Input: N qubits in |0...0⟩ state.
    // Goal: create an equal superposition of all basis vectors from |0...0⟩ to |1...1⟩
    // (i.e. state (|0...0⟩ + ... + |1...1⟩) / sqrt(2^N) ).
    operation AllBasisVectorsSuperposition (qs : Qubit[]) : Unit {
        // Apply H to all qubits
        for q in qs {
            H(q);
        }
    }


    // Task 1.10. Superposition of all even or all odd numbers
    // Inputs:
    //      1) N qubits in |0...0⟩ state.
    //      2) A boolean isEven.
    // Goal: create a superposition of all even numbers on N qubits if isEven is true,
    //       or a superposition of all odd numbers on N qubits if isEven is false.
    // 
    // A basis state encodes an integer number using big-endian binary notation: 
    // state |01⟩ corresponds to the integer 1, and state |10⟩ - to the integer 2. 
    //
    // Example: for N = 2 and isEven = true the required state is (|00⟩ + |10⟩) / sqrt(2), 
    //      and for N = 2 and isEven = false - (|01⟩ + |11⟩) / sqrt(2).
    operation EvenOddNumbersSuperposition (qs : Qubit[], isEven : Bool) : Unit {
        let N = Length(qs);
    
        // Even
        for i in 0 .. N-2 {
            H(qs[i]);
        }
        // Odd
        if (not isEven) {
            X(qs[N-1]);
        }
    }
    

    // Task 1.11. Superposition of |0...0⟩ and given bit string
    // Inputs:
    //      1) N qubits in |0...0⟩ state
    //      2) bit string represented as Bool[]
    // Goal: create an equal superposition of |0...0⟩ and basis state given by the bit string.
    // Bit values false and true correspond to |0⟩ and |1⟩ states.
    // You are guaranteed that the qubit array and the bit string have the same length.
    // You are guaranteed that the first bit of the bit string is true.
    // Example: for bit string = [true, false] the qubit state required is (|00⟩ + |10⟩) / sqrt(2).
    operation ZeroAndBitstringSuperposition (qs : Qubit[], bits : Bool[]) : Unit {
        // The following lines enforce the constraints on the input that you are given.
        // You don't need to modify them. Feel free to remove them, this won't cause your code to fail.
        Fact(Length(bits) == Length(qs), "Arrays should have the same length");
        Fact(Head(bits), "First bit of the input bit string should be set to true");
        
        //
        H(qs[0]);

        for i in 1 .. Length(qs) - 1 {
            if (bits[i]) {
                CNOT(qs[0], qs[i]);
            }
        }
    }


    // Task 1.12. Superposition of two bit strings
    // Inputs:
    //      1) N qubits in |0...0⟩ state
    //      2) two bit string represented as Bool[]s
    // Goal: create an equal superposition of two basis states given by the bit strings.
    //
    // Bit values false and true correspond to |0⟩ and |1⟩ states.
    // Example: for bit strings [false, true, false] and [false, false, true]
    // the qubit state required is (|010⟩ + |001⟩) / sqrt(2).
    // You are guaranteed that the both bit strings have the same length as the qubit array,
    // and that the bit strings will differ in at least one bit.
    function FindFirstDiff (bits1 : Bool[], bits2 : Bool[]) : Int {
        mutable firstDiff = -1;
        for i in 0 .. (Length(bits1) - 1) {
            if (bits1[i] != bits2[i] and firstDiff == -1) {
                set firstDiff = i;
            }
        }
        return firstDiff;
    }
    operation TwoBitstringSuperposition (qs : Qubit[], bits1 : Bool[], bits2 : Bool[]) : Unit {
        // find the index of the first bit at which the bit strings are different
        let firstDiff = FindFirstDiff (bits1, bits2);
            
        // Hadamard corresponding qubit to create superposition
        H(qs[firstDiff]);
            
        // iterate through the bit strings again setting the final state of qubits
        for i in 0 .. (Length(qs) - 1) {
            if (bits1[i] == bits2[i]) {
                // if two bits are the same apply X or nothing
                if (bits1[i]) {
                    X(qs[i]);
                }
            } else {
                // if two bits are different, set their difference using CNOT
                if (i > firstDiff) {
                    CNOT(qs[firstDiff], qs[i]);
                    if (bits1[i] != bits1[firstDiff]) {
                        X(qs[i]);
                    }
                }
            }
        }
    }


    // Task 1.13*. Superposition of four bit strings
    // Inputs:
    //      1) N qubits in |0...0⟩ state
    //      2) four bit string represented as Bool[][] bits
    //         bits is an array of size 4 x N array which describes the bit strings as follows:
    //         bits[i] describes the i-th bit string and has N elements.
    //         All four bit strings will be distinct.
    //
    // Goal: create an equal superposition of the four basis states given by the bit strings.
    //
    // Example: for N = 3 and bits = [[false, true, false], [true, false, false], [false, false, true], [true, true, false]]
    //          the state you need to prepare is (|010⟩ + |100⟩ + |001⟩ + |110⟩) / 2.
    operation FourBitstringSuperposition (qs : Qubit[], bits : Bool[][]) : Unit {
        // Hint: remember that you can allocate extra qubits.

        let N = Length(qs);
        
        use anc = Qubit[2] {
            // Put two ancillas into equal superposition of 2-qubit basis states |+>
            ApplyToEachA(H, anc);
            
            // If bits[i][j] is true, we preform (ControlledOnInt(i, X))(anc, qs[j])
            // Set up the right pattern on the main qubits with control on ancillas
            for i in 0 .. 3 {
                for j in 0 .. N - 1 {
                    if ((bits[i])[j]) {
                        (ControlledOnInt(i, X))(anc, qs[j]);
                    }
                }
            }
            
            // Then we have to iterate through each bitstring. 
            // If we are on bits[i] and i is odd, we do (ControlledOnBitString(bits[i], X))(qs, anc[0])
            // If it's even, we do (ControlledOnBitString(bits[i], X))(qs, anc[1])
            // Uncompute the ancillas, using patterns on main qubits as control
            for (i in 0 .. 3) {
                if (i % 2 == 1) {
                    (ControlledOnBitString(bits[i], X))(qs, anc[0]);
                }
                if (i / 2 == 1) {
                    (ControlledOnBitString(bits[i], X))(qs, anc[1]);
                }
            }
            
        }
    }

    // Solution 2
//     operation FourBitstringSuperposition (qs : Qubit[], bits : Bool[][]) : Unit {
//     use anc = Qubit[2];
//     // Put two ancillas into equal superposition of 2-qubit basis states
//     ApplyToEachA(H, anc);

//     // Set up the right pattern on the main qubits with control on ancillas
//     for i in 0 .. 3 {
//         for j in 0 .. Length(qs) - 1 {
//             if (bits[i][j]) {
//                 (ControlledOnInt(i, X))(anc, qs[j]);
//             }
//         }
//     }

//     // Uncompute the ancillas, using patterns on main qubits as control
//     for i in 0 .. 3 {
//         if (i % 2 == 1) {
//             (ControlledOnBitString(bits[i], X))(qs, anc[0]);
//         }
//         if (i / 2 == 1) {
//             (ControlledOnBitString(bits[i], X))(qs, anc[1]);
//         }
//     }
// }




    // Task 1.14. Superposition of all bit strings of the given parity
    // Inputs:
    //      1) N qubits in |0..0⟩ state (stored in an array of length N).
    //      2) An int "parity".
    // Goal: change the state to an equal superposition of all basis states that have
    //       an even number of 1s in them if "parity" = 0, or
    //       an odd number of 1s in them if "parity" = 1.
    // Example: for N = 2, the goal state would be (|00⟩ + |11⟩) / sqrt(2) if "parity" = 0,
    //       and (|01⟩ + |10⟩) / sqrt(2) if "parity" = 1.
    operation AllStatesWithParitySuperposition (qs : Qubit[], parity : Int) : Unit is Adj+Ctl {
        // Hint: remember that you can call the solution recursively.
        //       You are allowed to modify the signature of the method to include adjoint and/or controlled variants.
        
        // base of recursion: if N = 1, set the qubit to parity
        let N = Length(qs);
        if (N == 1) {
            if (parity == 1) {
                X(qs[0]);
            }
        } else {
            // split the first qubit into 0 and 1 (with equal amplitudes!)
            H(qs[0]);
            // prep 0 ⊗ state with the same parity and 1 ⊗ state with the opposite parity
            (ControlledOnInt(0, AllStatesWithParitySuperposition))(qs[0 .. 0], (qs[1 ...], parity));
            (ControlledOnInt(1, AllStatesWithParitySuperposition))(qs[0 .. 0], (qs[1 ...], 1 - parity));
        }
    }


    //////////////////////////////////////////////////////////////////
    // Part II. Arbitrary Rotations
    //////////////////////////////////////////////////////////////////

    // Task 2.1. Unequal superposition
    // Inputs:
    //      1) a qubit in the |0⟩ state.
    //      2) angle alpha, in radians, represented as Double.
    // Goal: prepare a cos(alpha) * |0⟩ + sin(alpha) * |1⟩ state on this qubit.
    operation UnequalSuperposition (q : Qubit, alpha : Double) : Unit {
        // Hint: Experiment with rotation gates from the Microsoft.Quantum.Intrinsic namespace.
        // Note that all rotation operators rotate the state by _half_ of its angle argument.
        Ry(2.0 * alpha, q);
    }


    // Task 2.2. 1/sqrt(2)|00⟩ + 1/2|01⟩ + 1/2|10⟩ state
    // Input: two qubits in |00⟩ state (stored in an array of length 2).
    // Goal: change the state to 1/sqrt(2)|00⟩ + 1/2|10⟩ + 1/2|11⟩.
    operation ControlledRotation (qs : Qubit[]) : Unit {
        // Put first qubit in alpha|0> + beta|1>
        H(qs[0]);
        
        //
        Controlled H ([qs[0]], qs[1]);    }


    // Task 2.3*. |00⟩ + |01⟩ + |10⟩ state
    // Input: 2 qubits in |00⟩ state (stored in an array of length 2).
    // Goal: change the state to (|00⟩ + |01⟩ + |10⟩) / sqrt(3).
    operation ThreeStates_TwoQubits (qs : Qubit[]) : Unit {
        // Rotate first qubit to (sqrt(2) |0⟩ + |1⟩) / sqrt(3) (task 1.4 from BasicGates kata)
        let theta = ArcSin(1.0 / Sqrt(3.0));
        Ry(2.0 * theta, qs[0]);

        // Split the state sqrt(2) |0⟩ ⊗ |0⟩ into |00⟩ + |01⟩
        (ControlledOnInt(0, H))([qs[0]], qs[1]); 
    }


    // Task 2.4*. (|00⟩ + ω |01⟩ + ω² |10⟩) / sqrt(3)
    // Input: two qubits in |00⟩ state (stored in an array of length 2).
    // Goal: change the state to (|00⟩ + ω |01⟩ + ω² |10⟩) / sqrt(3) where ω is exp(2πi/3).
    operation ThreeStates_TwoQubits_Phases (qs : Qubit[]) : Unit {
        // Prepare (1/sqrt(3))(|00> + |01> + |10>)
        ThreeStates_TwoQubits(qs);
        
        //  To get to the final state, we need to add the relative phases to 
        // both  |01> and |10> basis states without changing the |00> state.
        R1(2.0 * PI() / 3.0, qs[1]);
        R1(4.0 * PI() / 3.0, qs[0]);
    }


    // Task 2.5*. Hardy State
    // Input: 2 qubits in |00⟩ state.
    // Goal: create the state (3|00⟩ + |01⟩ + |10⟩ + |11⟩) / sqrt(12) on these qubits.
    operation Hardy_State (qs : Qubit[]) : Unit {
        // https://quantumcomputing.stackexchange.com/a/6844/2879
    
        Ry(2.0 * ArcCos(Sqrt(10.0 / 12.0)), qs[0]);
        (ControlledOnInt(0, Ry))([qs[0]], (2.0 * ArcCos(Sqrt(9.0 / 10.0)), qs[1]));
        (ControlledOnInt(1, Ry))([qs[0]], (2.0 * PI()/4.0, qs[1]));
        // In this special case the second Controlled Ry is equivalent to a Controlled H:
        // Controlled H([qs[0]], qs[1]);    
        }


    // Task 2.6*. W state on 2ᵏ qubits
    // Input: N = 2ᵏ qubits in |0...0⟩ state.
    // Goal: create a W state (https://en.wikipedia.org/wiki/W_state) on these qubits.
    // W state is an equal superposition of all basis states on N qubits of Hamming weight 1.
    // Example: for N = 4, W state is (|1000⟩ + |0100⟩ + |0010⟩ + |0001⟩) / 2.
    operation WState_PowerOfTwo (qs : Qubit[]) : Unit is Adj+Ctl {
        let N = Length(qs);

        if (N == 1) {
            // base of recursion: |1⟩
            X(qs[0]);
        } else {
            let K = N / 2;
            use anc = Qubit();
            H(anc);

            (ControlledOnInt(0, WState_PowerOfTwo))([anc], qs[0 .. K - 1]);
            (ControlledOnInt(1, WState_PowerOfTwo))([anc], qs[K .. N - 1]);

            for i in K .. N - 1 {
                CNOT(qs[i], anc);
            }
        }
    }


    // Task 2.7**. W state on arbitrary number of qubits
    // Input: N qubits in |0...0⟩ state (N is not necessarily a power of 2).
    // Goal: create a W state (https://en.wikipedia.org/wiki/W_state) on these qubits.
    // W state is an equal superposition of all basis states on N qubits of Hamming weight 1.
    // Example: for N = 3, W state is (|100⟩ + |010⟩ + |001⟩) / sqrt(3).
    operation WState_Arbitrary (qs : Qubit[]) : Unit {
        let N = Length(qs);
        Ry(2.0 * ArcSin(Sqrt(1.0/IntAsDouble(N))), qs[0]);
        for i in 1 .. N - 1 {
            (ControlledOnInt(0, Ry(2.0 * ArcSin(Sqrt(1.0/IntAsDouble(N - i))), _)))(qs[0 .. i-1], qs[i]);
        }
    }
}
