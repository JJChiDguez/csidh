# Stronger and Faster Side-Channel Protections for CSIDH

This is joined work with:

	* Daniel Cervantes-Vázquez,
	* Mathilde Chenu de La Morinerie, 
	* Jesús-Javier Chi-Domínguez, 
	* Luca De Feo,
	* Francisco Rodríguez-Henríquez, and
	* Benjamin Smith.

	This work uses the fact that the action used in the CSIDH protocol 
	can be computed with
			1. shortest differential addition chains,
			2. only Edwards curves (no hybrid between Montgomery and Edwards),
			3. A projective elligator with u randomly chosen from {2, (p-1)/2}, and
			4. the only two ``if statements'' used are for asking if a ``public'' 
			point is the infinity point.
	The implemented methods are Simba with one and two torsion points,
	which require dummy operations.
	In addition, a dummy-free CSIDH approach was implemented where each 
	secret exponent e_i corresponding to l_i is randomly chosen from
			1. the set of all odd integer number in {-b_i, ..., b_i} if b_i is odd, or
			2. the set of all even integer number in {-b_i, ..., b_i} if b_i is even.
	Here, b_i is the number of degree-(l_i) isogenies to be computed.

	[Notes]
	[1]
		The Simba method implemented is based on the work of Michael Meyer, 
		Fabio Campos, Steffen Reith: 
		"On Lions and Elligators: An efficient constant-time implementation of CSIDH". 
		IACR Cryptology ePrint Archive 2018: 1198 (2018).
	[2]
		This implementation re-use the field arithmetic implemented in the
		original CSIDH paper by Wouter Castryck, Tanja Lange, Chloe Martindale, 
		Lorenz Panny, Joost Renes: 
		"CSIDH: An Efficient Post-Quantum Commutative Group Action". 
		ASIACRYPT (3) 2018: 395-427.
		(its eprint version can be download at https://eprint.iacr.org/2018/383)

# ------------------------------------------------------------------------
This C code implementation was perfomed by:

	* Daniel Cervantes-Vázquez <dcervantes@computacion.cs.cinvestav.mx>,
	* Jesús-Javier Chi-Domínguez <jjchi@computacion.cs.cinvestav.mx, chidoys@gmail.com>, and
	* Francisco Rodríguez-Henríquez <francisco@cs.cinvestav.mx>.

# ------------------------------------------------------------------------
# C code
To compile the files you can do the following:

First, you can use any version of gcc compiler (just set the variable CC as 
an input of the Makefile [variable CC is optional, gcc is set by default]).

# Testing a CSIDH protocol (key exchange protocol)
[Compilation]

	(Using dummy operations and one torsion point)
		make csidh BITLENGTH_OF_P=512 TYPE=WITHDUMMY_1
	(Using dummy operations and two torsion points)
		make csidh BITLENGTH_OF_P=512 TYPE=WITHDUMMY_2
	(Dummy-free approach and using two torsion points)
		make csidh BITLENGTH_OF_P=512 TYPE=DUMMYFREE
		

[Execution]

		./bin/csidh


# Running-time: number of field operations
[Compilation]

	(Using dummy operations and one torsion point)
		make action_cost BITLENGTH_OF_P=512 TYPE=WITHDUMMY_1
	(Using dummy operations and two torsion points)
		make action_cost BITLENGTH_OF_P=512 TYPE=WITHDUMMY_2
	(Dummy-free approach and using two torsion points)
		make action_cost BITLENGTH_OF_P=512 TYPE=DUMMYFREE

[Execution]

		./bin/action_cost

# Running-time: number of clock cycles
[Compilation]

	(Using dummy operations and one torsion point)
		make action_timing BITLENGTH_OF_P=512 TYPE=WITHDUMMY_1
	(Using dummy operations and two torsion points)
		make action_timing BITLENGTH_OF_P=512 TYPE=WITHDUMMY_2
	(Dummy-free approach and using two torsion points)
		make action_timing BITLENGTH_OF_P=512 TYPE=DUMMYFREE

[Execution]

		./bin/action_timing

# Clean

	make clean
