 (* This defines the accessors for the parameter set of the dice model. It is stored as a list
    {type1Prob, type2Prob, faceProbs1, faceProbs2}. The accessors make the code more readable. *)
 type1Prob[parameters_] := parameters[[1]];
 type2Prob[parameters_] := parameters[[2]];
 faceProbs1[parameters_] := parameters[[3]];
 faceProbs2[parameters_] := parameters[[4]];
 
 (* diceEM just sets up for the actual EM algorithm by counting the frequencies of the faces
    in each trial and providing initial random or arbitrary values to diceEMIterator.*)
 (* A trial is the result of drawing a die at random from the bag and rolling it n times. *)
 diceEM[trials_, maxIterations_, accuracy_, randomSeed_:314, debug_:False]:=
	Module[{numFaces, binCountsList, initialFaceProbs1, initialFaceProbs2, nonNormalized1, nonNormalized2},
		(*SeedRandom[] initializes the random number generator so you can get deterministic
		  results for testing purposes.*)
		SeedRandom[randomSeed];
		(* For strictly defined EM, which uses maximum likelihood estimations, any die faces
		   that do not occur anywhere in the input can be treated as though they don't exist.*)
		numFaces = Max[trials];
		(*initialize binCountsList here, by mapping the builtin function BinCounts over the trials.
		  Each entry in the resulting binCountsList should be a list of the face frequencies from one trial.
		  E.g., if there are four faces on the dice, one entry might be {3, 0, 1, 1}, indicating 
		  that on this trial face 1 was rolled 3 times, face 2 was not rolled at all, etc.  *)
		  binCountsList=Table[BinCounts[trials[[i]], {1, numFaces+1}], {i, 1, Length[trials]}];
	    (* Initialize initialFaceProbs1 and initialFaceProbs2 here. Use RandomReal and make it so
	       so that all entries are within a factor of two of one another. The built in functions
	       Normalize and Total may also be useful here. *)
	    	nonNormalized1=Table[RandomReal[{1/numFaces, 2/numFaces}], {i, 1, numFaces}];
	    	initialFaceProbs1=Table[nonNormalized1[[i]]/Total[nonNormalized1], {i, 1, numFaces}];
	       
	   	   nonNormalized2=Table[RandomReal[{1/numFaces, 2/numFaces}], {i, 1, numFaces}];
	       initialFaceProbs2=Table[nonNormalized2[[i]]/Total[nonNormalized2], {i, 1, numFaces}];
	       
	diceEMIterator[binCountsList, 
		           numFaces, 
		           {0.45, 0.55, initialFaceProbs1, initialFaceProbs2}, 
		           maxIterations, 
		           accuracy]]
(* diceEMIterator implements the outer loop of the EM algorithm.
   It calls updateProbs on each iteration. *)		           
diceEMIterator[binCountsList_, numFaces_, initParams_, maxIterations_, accuracy_, debug_:False]:=
	Module[{oldParamEstimates, newParamEstimates},
		(* Initialize the local variables. *)
		oldParamEstimates = initParams;
		(* Loop here until either maxIterations has been reached or the accuracy goal has been
		   met. The accuracy goal is met when the sum, over all estimated parameters, of the absolute 
		   values of the changes from one iteration to the next is less than accuracy. *)
		(*Finally, if termination conditions have not been met, set old values to be the same 
		   as the new values.
		   *)
	     Do[
            newParamEstimates = updateProbs[binCountsList,oldParamEstimates];
            If[Total[Abs[newParamEstimates-oldParamEstimates],2]<=accuracy,
                Break[],
                oldParamEstimates = newParamEstimates
            ];,
            maxIterations
        ];
		(*At the end, return the estimated parameters with the less likely die first.*)
		If[type1Prob[newParamEstimates] <= type2Prob[newParamEstimates],
		   newParamEstimates,
		   {type2Prob[newParamEstimates], 
		    type1Prob[newParamEstimates], 
		    faceProbs2[newParamEstimates], 
		    faceProbs1[newParamEstimates]}]
	]
   
(* updateProbs does the actual EM calculations. *)
updateProbs[binCountsList_, oldParamEstimates_, debug_:False] :=
	Module[{posteriors,
		    (* type1Count and type2Count are the expected number of times a type1 or type2
		       die was drawn.*) 
		    type1Count, type2Count,
		    (* faceCounts1 is the expected number of times each face was rolled on a die 
		       of type 1.Likewise for faceCounts2.*) 
		    faceCounts1, faceCounts2, newFaceProbs1, newFaceProbs2},
		(*Create list of posterior probabilities of a Type1 die having been rolled on each draw 
		   by calling dicePosteriors. *)
		posteriors=Table[dicePosterior[binCountsList[[i]], type1Prob[oldParamEstimates], type2Prob[oldParamEstimates], faceProbs1[oldParamEstimates], faceProbs2[oldParamEstimates]],
			 	   {i,1,Length[binCountsList]}];
		(* Now use the posteriors to calculate EXPECTED number of times each die type was drawn. *) 
		type1Count=Sum[posteriors[[i]], {i, 1, Length[posteriors]}];
		type2Count=Sum[1-posteriors[[i]], {i, 1, Length[posteriors]}];
		(* Now use the posteriors to calculate EXPECTED number of times each face was rolled
		   on each die typep. *) 
		
		faceCounts1=Table[Sum[binCountsList[[i]][[j]]*posteriors[[i]], {i, 1, Length[binCountsList]}], {j, 1, Length[faceProbs1[oldParamEstimates]]}];
		faceCounts2=Table[Sum[binCountsList[[i]][[j]]*(1-posteriors[[i]]), {i, 1, Length[binCountsList]}], {j, 1, Length[faceProbs2[oldParamEstimates]]}];
		(* Finally, use these counts to compute maximum likelihood estimates for the parameters and 
		   return these estimates in a list: {newType1Prob, newType2Prob, newFaceProbs1, newFaceProbs2} *)
		newFaceProbs1=Table[faceCounts1[[i]]/Total[faceCounts1], {i, 1, Length[faceCounts1]}];
		newFaceProbs2=Table[faceCounts2[[i]]/Total[faceCounts2], {i, 1, Length[faceCounts2]}];
		
		{type1Count/(type2Count+type1Count), type2Count/(type2Count+type1Count), newFaceProbs1, newFaceProbs2}
	] 

(* diceSample creates a sample of rolls given the probabilities of dice 1, dice 2, the probabilities of their faces,
	 the rolls per draw and draws. *)
   diceSample[numType1_, numType2_, type1_, type2_, draws_, rollsPerDraw_] :=
 Module[{result=List[]},
 	(*Do draws times the following:*)
 	Do[
 	(*Append to our result the following:*)
 		AppendTo[result,
 	(*A random variate (which will be repeated rollsPerDraw times)...*)
 			RandomVariate[
 	(*...with an empirical distribution of either type2 or type1 which will depend on...*)
 				EmpiricalDistribution[
 	(*...weather for our particular role dice 1 or dice 2 is more probable (based on a bernoulli distribution), as written
 	if the distribution gives a variate of 1, then type 2 is more likely, else it is type 1.*)
 					If[RandomVariate[BernoulliDistribution[numType2/(numType2+numType1)]]==1,
 	(*So if it is type 2 or empirical distribution will follow the distribution of type2 with respect to the number of
 	sides it has (range of the length of the list)*)
 					type2 -> Range[Length[type2]],
 	(*We will follow the above steps if type 1 is more likely*)
 					type1 -> Range[Length[type1]]]],
 			rollsPerDraw]], 
 	draws];
 result]
 
 (*The following 4 methods (dicePosterior, probabilityGiven1, probabilityGiven2, NumOfSides) are for posterior probability*)
 	(*The posterior probability of type 1 is:*)
dicePosterior[binCounts_, type1Prior_, type2Prior_, faceProbs1_, faceProbs2_] := 
Module[{},
	(*The probability of our observations given dice 1, times the prior probability of dice 1*)
	(probabilityGiven1[binCounts, faceProbs1]*type1Prior)
	(*Divided by the total probability of our observations, which is calculated by exhaustive conditionalization on
	our probabilities given both die and the priors.*)
	/(type1Prior*probabilityGiven1[binCounts, faceProbs1]+type2Prior*probabilityGiven2[binCounts, faceProbs2])
	]

probabilityGiven1[binCounts_, faceProbs1_]:= 
Module[{binSize=NumOfSides[binCounts], result=1},
	(*The probability of our results given dice 1 are the probability of the face at each "index" to the number of times
	this face was seen.*)
Do[If[faceProbs1[[i]]==0 && binCounts[[i]]==0, result=result, result=result*(faceProbs1[[i]]^binCounts[[i]])], {i, binSize}]; result]
	
probabilityGiven2[binCounts_, faceProbs2_]:= 
Module[{binSize=NumOfSides[binCounts], result=1},
Do[If[faceProbs2[[i]]==0 && binCounts[[i]]==0, result=result, result=result*(faceProbs2[[i]]^binCounts[[i]])], {i, binSize}]; result]

NumOfSides[faceProbs1_]:= Length[faceProbs1]

(*myRound is a hack to get around a problem with rounding numbers in \
Mathematica.*)

myRound[x_, n_] :=
  N[IntegerPart[Round[x,10^-n]*10^n] / 10^n];