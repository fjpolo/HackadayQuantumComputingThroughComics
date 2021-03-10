// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT license.

//////////////////////////////////////////////////////////////////
// This file is a back end for the tasks in Deutsch-Jozsa algorithm tutorial.
// We strongly recommend to use the Notebook version of the tutorial
// to enjoy the full experience.
//////////////////////////////////////////////////////////////////

namespace Quantum.Kata.DeutschJozsaAlgorithm {
    
    open Microsoft.Quantum.Diagnostics;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Canon;
    
    
    //////////////////////////////////////////////////////////////////
    // Part I. Classical algorithm
    //////////////////////////////////////////////////////////////////

    // Exercise 1.
    function Function_MostSignificantBit (x : Int, N : Int) : Int {
        return x >>> (N-1);
    }


    // Exercise 2. 
    function IsFunctionConstant_Classical (N : Int, f : (Int -> Int)) : Bool {
    // Store first value
    let firstVal = f(0);
    // Try all the following inputs to see if any of the values differ
    for input in 1..2 <<< (N-1){
        let nextVal = f(input);
        if nextVal != firstVal{
            // Two different values -> Balanced
            return false;
        } // END if
    } // END for
    // Got 2ᴺ⁻¹ + 1 copies of the same value -> Constant
    return true;
} // END IsFunctionConstant_Classical()


    //////////////////////////////////////////////////////////////////
    // Part II. Quantum oracles
    //////////////////////////////////////////////////////////////////
    
    // Exercise 3.
    operation PhaseOracle_MostSignificantBit (x : Qubit[]) : Unit {
        Z(x[0]);
    }


    //////////////////////////////////////////////////////////////////
    // Part III. Quantum algorithm
    //////////////////////////////////////////////////////////////////
    
    // Exercise 4.
    operation DeutschJozsaAlgorithm (N : Int, oracle : (Qubit[] => Unit)) : Bool {
        // Create a boolean variable for storing the return value.
        // You'll need to update it later, so it has to be declared as mutable.
        mutable isConstant = true;

        // Allocate an array of N qubits for the input register x.
        use x = Qubit[N];
        // Newly allocated qubits start in the |0⟩ state.

        // The first step is to prepare the qubits in the required state before calling the oracle.
        // A qubit can be transformed from the |0⟩ state to the |+⟩ state by applying a Hadamard gate H.
        // To apply this to each qubit, you can use a for loop to iterate over all array elements
        // using the following syntax: for q in qs { ... }
        ApplyToEach(H, x);

        // Apply the oracle to the input register.
        // The syntax is the same as for applying any function or operation.
        oracle(x);

        // Apply a Hadamard gate to each qubit of the input register again.
        ApplyToEach(H, x);

        // Measure each qubit of the input register in the computational basis using the M operation.
        // You can use a for loop to iterate over the qubits of the array again.
        // Note that you can't return the answer in the middle of a loop,
        // you have to update the variable isConstant using the "set" keyword.
        for q in x {
            if (M(q) == One) {
                set isConstant = false;
            }
        }
        
        // Return the value of the boolean variable.
        return isConstant;
    }
}