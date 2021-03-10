// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT license.

namespace Quantum.Kata.Measurements {

    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Convert;
    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Measurement; 
    open Microsoft.Quantum.Arithmetic;


    //////////////////////////////////////////////////////////////////
    // Welcome!
    //////////////////////////////////////////////////////////////////

    // The "Measurements" quantum kata is a series of exercises designed
    // to get you familiar with programming in Q#.
    // It covers the following topics:
    //  - using single-qubit measurements,
    //  - discriminating orthogonal and nonorthogonal states.

    // Each task is wrapped in one operation preceded by the description of the task.
    // Each task (except tasks in which you have to write a test) has a unit test associated with it,
    // which initially fails. Your goal is to fill in the blank (marked with // ... comment)
    // with some Q# code to make the failing test pass.

    // The tasks are given in approximate order of increasing difficulty; harder ones are marked with asterisks.

    //////////////////////////////////////////////////////////////////
    // Part I. Discriminating Orthogonal States
    //////////////////////////////////////////////////////////////////

    // Task 1.1. |0⟩ or |1⟩ ?
    // Input: a qubit which is guaranteed to be in either the |0⟩ or the |1⟩ state.
    // Output: true if the qubit was in the |1⟩ state, or false if it was in the |0⟩ state.
    // The state of the qubit at the end of the operation does not matter.
    operation IsQubitOne (q : Qubit) : Bool {
        // The operation M will measure a qubit in the Z basis (|0⟩ and |1⟩ basis)
        // and return Zero if the observed state was |0⟩ or One if the state was |1⟩.
        // To answer the question, you need to perform the measurement and check whether the result
        // equals One - either directly or using library function IsResultOne.
        //
        // Replace the returned expression with (M(q) == One).
        // Then rebuild the project and rerun the tests - T101_IsQubitOne should now pass!
        return M(q) == One;
    }


    // Task 1.2. Set qubit to |0⟩ state
    // Input: a qubit in an arbitrary state.
    // Goal:  change the state of the qubit to |0⟩.
    operation InitializeQubit (q : Qubit) : Unit {
        // if q == |1> -> X(q) = |0>
        if(M(q) == One){
            X(q);
        }
    }


    // Task 1.3. |+⟩ or |-⟩ ?
    // Input: a qubit which is guaranteed to be in either the |+⟩ or the |-⟩ state
    //        (|+⟩ = (|0⟩ + |1⟩) / sqrt(2), |-⟩ = (|0⟩ - |1⟩) / sqrt(2)).
    // Output: true if the qubit was in the |+⟩ state, or false if it was in the |-⟩ state.
    // The state of the qubit at the end of the operation does not matter.
    operation IsQubitPlus (q : Qubit) : Bool {
        // H(|0>) = |+> ^ H(|+>) = |0>
        // H(|0>) = |+> ^ H(|+>) = |0>
        // H(|1>) = |-> ^ H(|->) = |1>
        H(q);
        //
        //if(M(q) == One){
        //    return false;
        //}
        //else{
        //    return true;
        //}
        return M(q) == Zero;
    }


    // Task 1.4. |A⟩ or |B⟩ ?
    // Inputs:
    //      1) angle α, in radians, represented as a Double
    //      2) a qubit which is guaranteed to be in either the |A⟩ or the |B⟩ state, where
    //         |A⟩ =   cos α |0⟩ + sin α |1⟩,
    //         |B⟩ = - sin α |0⟩ + cos α |1⟩.
    // Output: true if the qubit was in the |A⟩ state, or false if it was in the |B⟩ state.
    // The state of the qubit at the end of the operation does not matter.
    operation IsQubitA (alpha : Double, q : Qubit) : Bool {
        Ry(-2.0 * alpha, q);
        return M(q) == Zero;
    }


    // Task 1.5. |00⟩ or |11⟩ ?
    // Input: two qubits (stored in an array of length 2) which are guaranteed to be in either the |00⟩ or the |11⟩ state.
    // Output: 0 if the qubits were in the |00⟩ state,
    //         1 if they were in |11⟩ state.
    // The state of the qubits at the end of the operation does not matter.
    operation ZeroZeroOrOneOne (qs : Qubit[]) : Int {
        //  Since there are no possible mixed states, we only have to measure 
        // one of the qubits to determine the state of the other.
        if(M(qs[0]) == One){
            return 1;
        }
        else{
            return 0;
        }
    }


    // Task 1.6. Distinguish four basis states
    // Input: two qubits (stored in an array) which are guaranteed to be
    //        in one of the four basis states (|00⟩, |01⟩, |10⟩ or |11⟩).
    // Output: 0 if the qubits were in |00⟩ state,
    //         1 if they were in |01⟩ state,
    //         2 if they were in |10⟩ state,
    //         3 if they were in |11⟩ state.
    // In this task and the subsequent ones the order of qubit states
    // in task description matches the order of qubits in the array
    // (i.e., |10⟩ state corresponds to qs[0] in state |1⟩ and qs[1] in state |0⟩).
    // The state of the qubits at the end of the operation does not matter.
    operation BasisStateMeasurement (qs : Qubit[]) : Int {
        if(M(qs[0]) == Zero){
            if(M(qs[1]) == Zero){
                return 0;
            }
            else{
                return 1;
            }
        }
        else{
            if(M(qs[1]) == Zero){
                return 2;
            }
            else{
                return 3;
            }    
        }
    }


    // Task 1.7. Distinguish two basis states given by bit strings
    // Inputs:
    //      1) N qubits (stored in an array) which are guaranteed to be
    //         in one of the two basis states described by the given bit strings.
    //      2) two bit string represented as Bool[]s.
    // Output: 0 if the qubits were in the basis state described by the first bit string,
    //         1 if they were in the basis state described by the second bit string.
    // Bit values false and true correspond to |0⟩ and |1⟩ states.
    // The state of the qubits at the end of the operation does not matter.
    // You are guaranteed that both bit strings have the same length as the qubit array,
    // and that the bit strings differ in at least one bit.
    // You can use exactly one measurement.
    // Example: for bit strings [false, true, false] and [false, false, true]
    //          return 0 corresponds to state |010⟩, and return 1 corresponds to state |001⟩.
   function FindFirstDiff (bits1 : Bool[], bits2 : Bool[]) : Int {
        for i in 0 .. (Length(bits1) - 1) {
            if (bits1[i] != bits2[i]) {
                return i;
            }
        }
    return -1;
    }
    
    operation TwoBitstringsMeasurement (qs : Qubit[], bits1 : Bool[], bits2 : Bool[]) : Int {
        // First, we need to find where the bitstrings diverge.
        let firstDiff = FindFirstDiff(bits1, bits2);
        
        //  Then we can measure our qubit array at the FirstDiff index. We can turn the result 
        // of our measurement into a boolean by testing it's equality to One.
        //
        //  If our qubit measured |1⟩, our variable will be true. If we measured |0⟩, our variable 
        // will be false.
        //
        //  Now we can check if the bit we measured belongs to the state describe by the first or 
        // second bitstring.
        //
        //  If the boolean at bits1[firstDiff] is the same as our measured value, we know bits1[] 
        // describes our current state and we can return 0.
        
        // res = true if the first different qubit measures to be one
        let res = M(qs[firstDiff]) == One;
        
        // If the measurement aligns with the value in the first bitstring, return 0
        // Otherwise it's the second state, return false.
        return res == bits1[firstDiff] ? 0 | 1;
    }


    // Task 1.8. Distinguish two superposition states given by two arrays of bit strings - 1 measurement
    // Inputs:
    //      1) N qubits which are guaranteed to be
    //         in one of the two superposition states described by the given arrays of bit strings.
    //      2) two arrays of bit strings represented as Bool[][]s.
    //         The arrays have dimensions M₁ ⨯ N and M₂ ⨯ N respectively, where N is the number of
    //         qubits and M₁ and M₂ are the numbers of bit strings in each array. Note that in general M₁ ≠ M₂.
    //         An array of bit strings [b₁, ..., bₘ] defines a state that is
    //         an equal superposition of all basis states defined by bit strings b₁, ..., bₘ.
    //         For example, an array of bit strings [[false, true, false], [false, true, true]]
    //         defines a superposition state (|010⟩ + |011⟩) / sqrt(2).
    //         
    // Output: 0 if qubits were in the superposition state described by the first array,
    //         1 if they were in the superposition state described by the second array.
    // The state of the qubits at the end of the operation does not matter.
    //
    // You are allowed to use exactly one measurement.
    // You are guaranteed that there exists an index of a qubit Q for which 
    //  - all the bit strings in the first array have the same value in this position (all bits1[j][Q] are the same),
    //  - all the bit strings in the second array have the same value in this position (all bits2[j][Q] are the same),
    //  - these values are different for the first and the second arrays.
    // 
    // Example: for arrays [[false, true, false], [false, true, true]] and [[true, false, true], [false, false, true]]
    //          return 0 corresponds to state (|010⟩ + |011⟩) / sqrt(2), 
    //          return 1 corresponds to state (|101⟩ + |001⟩) / sqrt(2),
    //          and you can distinguish these states perfectly by measuring the second qubit.
    
    function FindFirstSuperpositionDiff ( bits1:Bool[][], bits2:Bool[][], nQubits:Int ) :Int {
        //
        //  Properties:
        //    1- The value of all arrays in bits1 at the index Q is either true or false, 
        //      and the value of all arrays in bits2 at the index Q is either true or false. 
        //      If this is not the case, you cannot be sure that measuring the corresponding 
        //      qubit will always return the same result.
        //    2- This value is different for bits1 and bits2.
        //
        //  We will iterate over all qubit indices, and for each of them we'll calculate the 
        // number of 1s in that position in bits1 and bits2. 
        //    1. The first condition means that this count should equal 0 (if all bit strings 
        //      have 0 bit in this position) or the length of the array of bit strings (if all 
        //      bit strings have 1 bit in this position). 
        //    2. The second condition means that this count is different for bits1 and bits2, 
        //      i.e., one of the counts should equal 0 and another one - the length of the 
        //      corresponding array of bit strings.
        //
        
        // Iterate through qubit indices
        for i in 0 .. (nQubits-1){
            // count the number of 1s in i-th position in bit strings of both arrays
            mutable val1 = 0;
            mutable val2 = 0;
            
            // Iterate through bits1
            for j in 0 .. (Length(bits1)-1) {
                if bits1[j][i]{
                    set val1 += 1;
                }        
            } // END for bits1
            
            // Iterate through bits2
            for k in 0 .. (Length(bits2)-1) {
                if bits2[k][i]{
                    set val2 += 1;
                } 
            } // END for bits2
            
            // Now compare val1 and val2 to get the position
            if ( ( val1 == Length(bits1) and val2 == 0 ) or 
                ( val1 == 0 and val2 == Length(bits2) )
            ) {
                return i;
            }
        } // END for
        
        // No match
        return -1;    
    } // END FindFirstSuperpositionDiff()
    operation SuperpositionOneMeasurement (qs : Qubit[], bits1 : Bool[][], bits2 : Bool[][]) : Int {
        //  We are looking for the index Q where the two bit strings differ. 
        // Let's define a function FindFirstSuperpositionDiff()
        let diff = FindFirstSuperpositionDiff(bits1, bits2, Length(qs));
            
        // Given the index we just found, we measure the qubit on that position.
        let res = ResultAsBool(M(qs[diff]));
        
        // Compare and return
        if (res == bits1[0][diff]) {
            return 0;
        }
        else {
            return 1;
        }
    }


    // Task 1.9. Distinguish two superposition states given by two arrays of bit strings
    // Inputs:
    //      1) N qubits which are guaranteed to be
    //         in one of the two superposition states described by the given arrays of bit strings.
    //      2) two arrays of bit strings represented as Bool[][]s.
    //         The arrays describe the superposition states in the same way as in the previous task,
    //         i.e., they have dimensions M₁ ⨯ N and M₂ ⨯ N respectively, N being the number of qubits.
    //
    // Output: 0 if qubits were in the superposition state described by the first array,
    //         1 if they were in the superposition state described by the second array.
    // The state of the qubits at the end of the operation does not matter.
    //
    // You can use as many measurements as you wish.
    // The only constraint on the bit strings is that all bit strings in the two arrays are distinct. 
    //
    // Example: for arrays [[false, true, false], [false, false, true]] and [[true, true, true], [false, true, true]]
    //          return 0 corresponds to state (|010⟩ + |001⟩) / sqrt(2), 
    //          return 1 corresponds to state (|111⟩ + |011⟩) / sqrt(2)
    operation SuperpositionMeasurement (qs : Qubit[], bits1 : Bool[][], bits2 : Bool[][]) : Int {
        //  When we measure all qubits of a certain superposition state, it collapses to one 
        // of the basis vectors that comprised the superposition. We can do exactly that and 
        // compare the resulting state to the given bit strings to see which array it belongs to.
        
        // Measure all qubits and store in array
        let measuredState = ResultArrayAsInt(MultiM(qs));
        
        // Iterate through bytes1
        for s in bits1 {
            // if measured state is the same as in bytes1, return 0 
            // meaning superposition state from first array
            if (BoolArrayAsInt(s) == measuredState) {
                return 0;
            }
        }
        return 1;
    }


    // Task 1.10. |0...0⟩ state or W state ?
    // Input: N qubits (stored in an array) which are guaranteed to be
    //        either in the |0...0⟩ state
    //        or in the W state (https://en.wikipedia.org/wiki/W_state).
    // Output: 0 if the qubits were in the |0...0⟩ state,
    //         1 if they were in the W state.
    // The state of the qubits at the end of the operation does not matter.
    operation AllZerosOrWState (qs : Qubit[]) : Int {
        //  Every possible measurement in a W state includes a single 1. If we 
        // find even a single one in our measured state, it can no longer be in 
        // state |0...0⟩ and thus it must be a W state.
        
        // Iterate through qubits
        for q in 0 .. (Length(qs)-1){
            if ResultAsBool(M(qs[q])) {
                return 1; // We found a one, it's W
            }
        }
        return 0; // it's |0...0⟩
    }


    // Task 1.11. GHZ state or W state ?
    // Input: N >= 2 qubits (stored in an array) which are guaranteed to be
    //        either in the GHZ state (https://en.wikipedia.org/wiki/Greenberger%E2%80%93Horne%E2%80%93Zeilinger_state)
    //        or in the W state (https://en.wikipedia.org/wiki/W_state).
    // Output: 0 if the qubits were in the GHZ state,
    //         1 if they were in the W state.
    // The state of the qubits at the end of the operation does not matter.
    operation GHZOrWState (qs : Qubit[]) : Int {
        //  If our qubits are in a GHZ state, we can expect them to be either 
        // all |0⟩ or all |1⟩. If we find two of our qubits measure as different 
        // values, we can assume they are in the W state.
        
        // Loop through qubits
        for q in 1 .. (Length(qs)-1){
            if ResultAsBool(M(qs[q-1])) != ResultAsBool(M(qs[q])){
                // W state
                return 1;
            }
        }
        // GHZ state
        return 0;
    }


    // Task 1.12. Distinguish four Bell states
    // Input: two qubits (stored in an array) which are guaranteed to be in one of the four Bell states:
    //         |Φ⁺⟩ = (|00⟩ + |11⟩) / sqrt(2)
    //         |Φ⁻⟩ = (|00⟩ - |11⟩) / sqrt(2)
    //         |Ψ⁺⟩ = (|01⟩ + |10⟩) / sqrt(2)
    //         |Ψ⁻⟩ = (|01⟩ - |10⟩) / sqrt(2)
    // Output: 0 if the qubits were in |Φ⁺⟩ state,
    //         1 if they were in |Φ⁻⟩ state,
    //         2 if they were in |Ψ⁺⟩ state,
    //         3 if they were in |Ψ⁻⟩ state.
    // The state of the qubits at the end of the operation does not matter.
    operation BellState (qs : Qubit[]) : Int {
        // Hint: you need to use 2-qubit gates to solve this task

        //  The first step is to apply the H() gate to both of our qubits.
        H(qs[0]);
        H(qs[1]);

        //  Then, we can reverse our Bell circuit by preforming a CNOT 
        // using the second bit as control, and then applying the 
        // H() gate to it.
        CNOT(qs[1], qs[0]);
        H(qs[1]);

        //  This will change our state to one of the following: 
        //        |00⟩, |01⟩, |10⟩, |11⟩.

        // We can then create two integer variables my measuring each 
        // of the qubits and testing their equality to Zero or One.
        let m1 = M(qs[0]) == Zero ? 0 | 1;
        let m2 = M(qs[1]) == Zero ? 0 | 1;

        // We can finally compute m2 * 2 + m1 and return it.
        return m2 * 2 + m1;
    }


    // Task 1.13. Distinguish four orthogonal 2-qubit states
    // Input: two qubits (stored in an array) which are guaranteed to be in one of the four orthogonal states:
    //         |S0⟩ = (|00⟩ + |01⟩ + |10⟩ + |11⟩) / 2
    //         |S1⟩ = (|00⟩ - |01⟩ + |10⟩ - |11⟩) / 2
    //         |S2⟩ = (|00⟩ + |01⟩ - |10⟩ - |11⟩) / 2
    //         |S3⟩ = (|00⟩ - |01⟩ - |10⟩ + |11⟩) / 2
    // Output: 0 if qubits were in |S0⟩ state,
    //         1 if they were in |S1⟩ state,
    //         2 if they were in |S2⟩ state,
    //         3 if they were in |S3⟩ state.
    // The state of the qubits at the end of the operation does not matter.
    operation TwoQubitState (qs : Qubit[]) : Int {
        //  We can apply the Hadamard gate to both our 
        // qubits to collapse them into one of the two qubit 
        // basis states |00⟩, |01⟩, |10⟩, or |11⟩.
        H(qs[0]);
        H(qs[1]);
        
        //  Then we can use the solution from task 1.6 to identify 
        // the correct basis state.
        return BasisStateMeasurement(qs);
    }


    // Task 1.14*. Distinguish four orthogonal 2-qubit states, part two
    // Input: two qubits (stored in an array) which are guaranteed to be in one of the four orthogonal states:
    //         |S0⟩ = ( |00⟩ - |01⟩ - |10⟩ - |11⟩) / 2
    //         |S1⟩ = (-|00⟩ + |01⟩ - |10⟩ - |11⟩) / 2
    //         |S2⟩ = (-|00⟩ - |01⟩ + |10⟩ - |11⟩) / 2
    //         |S3⟩ = (-|00⟩ - |01⟩ - |10⟩ + |11⟩) / 2
    // Output: 0 if qubits were in |S0⟩ state,
    //         1 if they were in |S1⟩ state,
    //         2 if they were in |S2⟩ state,
    //         3 if they were in |S3⟩ state.
    // The state of the qubits at the end of the operation does not matter.
    operation TwoQubitStatePartTwo (qs : Qubit[]) : Int {
        //  The trick to this task is applying the Hadamard gate to 
        // transform the four orthogonal states into the Bell states.
        H(qs[1]);

        // After this is achieved, we can implement code similar to 
        // task 1.10 to determine what the bell state is.
        CNOT(qs[0], qs[1]);
        H(qs[0]);
        let m1 = M(qs[0]) == One ? 0 | 1;
        let m2 = M(qs[1]) == One ? 0 | 1;
        return m2 * 2 + m1;
    }


    // Task 1.15**. Distinguish two orthogonal states on three qubits
    // Input: Three qubits (stored in an array) which are guaranteed to be in either one of the
    //        following two states:
    //        1/sqrt(3) ( |100⟩ + ω |010⟩ + ω² |001⟩ ),
    //        1/sqrt(3) ( |100⟩ + ω² |010⟩ + ω |001⟩ ).
    //        Here ω = exp(2π i/3) denotes a primitive 3rd root of unity.
    // Output: 0 if the qubits were in the first superposition,
    //         1 if they were in the second superposition.
    // The state of the qubits at the end of the operation does not matter.
    operation WState_Arbitrary (qs : Qubit[]) : Unit is Adj + Ctl {
        let N = Length(qs);

        if (N == 1) {
            // base case of recursion: |1⟩
            X(qs[0]);
        } else {
            // |W_N⟩ = |0⟩|W_(N-1)⟩ + |1⟩|0...0⟩
            // do a rotation on the first qubit to split it into |0⟩ and |1⟩ with proper weights
            // |0⟩ -> sqrt((N-1)/N) |0⟩ + 1/sqrt(N) |1⟩
            let theta = ArcSin(1.0 / Sqrt(IntAsDouble(N)));
            Ry(2.0 * theta, qs[0]);

            // do a zero-controlled W-state generation for qubits 1..N-1
            X(qs[0]);
            Controlled WState_Arbitrary(qs[0 .. 0], qs[1 .. N - 1]);
            X(qs[0]);
        }
    }
    operation ThreeQubitMeasurement (qs : Qubit[]) : Int {
        //  In order to solve this, we must map our first state to a 
        // W state, and our second state to the first.We can achieve 
        // this by applying a unitary operation of the form I2⊗R⊗R2.
        //
        // Now we can take the new W state and undo it into the state 
        // |000⟩, by utilizing the Adjoint if the operation that let us 
        // prepare the W state in the last kata.
        // 
        // Our second state will then be mapped to some state guaranteed 
        // to be perpendicular to |000⟩. In simple terms, the second state 
        // will never be measured to be |000⟩.
        //
        // Now we can measure our state, which is either |000⟩ (state 1) or 
        // something else (state 2). 
        
        //  Map the first state to 000 state and the second one to something orthogonal to it
        // (as described in reference solution)
        R1(-2.0 * PI() / 3.0, qs[1]);
        R1(-4.0 * PI() / 3.0, qs[2]);
        Adjoint WState_Arbitrary(qs);
        
        //  Measure all qubits: if all of them are 0, we have the first state,
        // if at least one of them is 1, we have the second state
        return MeasureInteger(LittleEndian(qs)) == 0 ? 0 | 1;
    }


    //////////////////////////////////////////////////////////////////
    // Part II*. Discriminating Nonorthogonal States
    //////////////////////////////////////////////////////////////////

    // The solutions for tasks in this section are validated using the following method.
    // The solution is called on N input states, each of which is picked randomly,
    // with all possible input states equally likely to be generated.
    // The accuracy of state discrimination is estimated as an average of
    // discrimination correctness over all input states.

    // Task 2.1*. |0⟩ or |+⟩ ?
    //           (quantum hypothesis testing or state discrimination with minimum error)
    // Input: a qubit which is guaranteed to be in either the |0⟩ or the |+⟩ state with equal probability.
    // Output: true if qubit was in the |0⟩ state, or false if it was in the |+⟩ state.
    // The state of the qubit at the end of the operation does not matter.
    // Note: in this task you have to get accuracy of at least 80%.
    operation IsQubitPlusOrZero (q : Qubit) : Bool {
        // ...
        return true;
    }


    // Task 2.2**. |0⟩, |+⟩ or inconclusive?
    //             (unambiguous state discrimination)
    // Input: a qubit which is guaranteed to be in either the |0⟩ or the |+⟩ state with equal probability.
    // Output: 0 if qubit was in the |0⟩ state,
    //         1 if it was in the |+⟩ state,
    //         -1 if you can't decide, i.e., an "inconclusive" result.
    // Your solution:
    //  - should never give 0 or 1 answer incorrectly (i.e., identify |0⟩ as 1 or |+⟩ as 0).
    //  - may give an inconclusive (-1) answer in at most 80% of the cases.
    //  - must correctly identify |0⟩ state as 0 in at least 10% of the cases.
    //  - must correctly identify |+⟩ state as 1 in at least 10% of the cases.
    //
    // The state of the qubit at the end of the operation does not matter.
    operation IsQubitPlusZeroOrInconclusiveSimpleUSD (q : Qubit) : Int {
        // ...
        return -2;
    }


    // Task 2.3**. Unambiguous state discrimination of 3 non-orthogonal states on one qubit
    //             (a.k.a. the Peres/Wootters game)
    // Input: a qubit which is guaranteed to be in one of the three states with equal probability:
    //        |A⟩ = 1/sqrt(2) (|0⟩ + |1⟩),
    //        |B⟩ = 1/sqrt(2) (|0⟩ + ω |1⟩),
    //        |C⟩ = 1/sqrt(2) (|0⟩ + ω² |1⟩),
    //          where ω = exp(2iπ/3) denotes a primitive, complex 3rd root of unity.
    // Output: 1 or 2 if the qubit was in the |A⟩ state,
    //         0 or 2 if the qubit was in the |B⟩ state,
    //         0 or 1 if the qubit was in the |C⟩ state.
    // The state of the qubit at the end of the operation does not matter.
    // You can use extra qubit(s) in your solution.
    // Note: in this task you have to succeed with probability 1, i.e., you are never allowed
    //       to give an incorrect answer.
    operation IsQubitNotInABC (q : Qubit) : Int {
        // ...
        return -1;
    }
}
