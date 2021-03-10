// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT license.

//////////////////////////////////////////////////////////////////
// This file is a back end for the tasks in the tutorial.
// We strongly recommend to use the Notebook version of the tutorial
// to enjoy the full experience.
//////////////////////////////////////////////////////////////////

namespace Quantum.Kata.SingleQubitSystemMeasurements {
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Math;
    
    // Exercise 2. Distinguish |0❭ and |1❭
    operation IsQubitZero (q : Qubit) : Bool {
        return M(q) == Zero;
    }

    // Exercise 3. Distinguish |+❭ and |-❭ using Measure operation
    operation IsQubitMinus (q : Qubit) : Bool {
        // https://www.quantum-inspire.com/kbase/qubit-basis-states/
        // |+> -> X -> 1
        // |-> -> -X -> -1
        return Measure([PauliX], [q]) == One;
    }

    // Exercise 5. Distinguish specific orthogonal states
    operation IsQubitPsiPlus (q : Qubit) : Bool {
        //  We can distinguish between the states |ψ±⟩ if we implement 
        // a measurement in the {|ψ+⟩, |ψ−⟩} basis. This can be done if 
        // we construct a unitary transformation which maps the |ψ+⟩ 
        // state to the |0⟩ state, and the |ψ−⟩ state to the |1⟩ state.
    
        // Map |ψ+> to |0> and |ψ-> to |1>
        // Apply Ry(−θ)
        Ry(-2.0 * ArcTan2(0.8, 0.6), q);
        // Measure
        return M(q) == Zero;
    }

    // Exercise 6. Distinguish states |A❭ and |B❭
    operation IsQubitA (alpha : Double, q : Qubit) : Bool {
        //  We can distinguish between the states |A⟩ and |B⟩ 
        // if we implement a measurement in the {|A⟩, |B⟩} basis.
        
        // Transform |A> to |0> and |B> to |1>
        // Apply Rx(−θ)
        Rx(-2.0 * alpha, q);
        // Measure
        return M(q) == Zero;
    }

    // Exercise 7. Measure state in {|A❭, |B❭} basis
    operation MeasureInABBasis (alpha : Double, q : Qubit) : Result {
        //  We first apply Rx(−θ) to the qubit. Next, we measure 
        // in the computational basis using the M operation. If 
        // the M operation returned Zero, we get measurement outcome  
        // A, and if it returned One, we get measurement outcome B. 
        // 
        //  After the measurement, we apply the inverse of the Rx(−θ) 
        // gate, which is the Rx(θ) gate. The final rotation ensures 
        // that the state of the qubit is in the state corresponding 
        // to the measurement outcome.
        Rx(-2.0 * alpha, q);
        let measurementResult = M(q);
        Rx(2.0 * alpha, q);
        return measurementResult;
    }
}

