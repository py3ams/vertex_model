function compile_mex_functions(compile_mex_functions_logical)

if compile_mex_functions_logical

	mex CalculateCellAreas.cpp
	mex CalculateCellConcentrations.cpp
	mex CalculateCentroid.cpp
	mex CalculateIngestionTerm.cpp
	mex CalculateSourceTerm.cpp
	mex CalculateTotalDpp.cpp
	mex CalculateTotalDppConst.cpp
	mex CreateBoundaryElement.cpp
	mex CreateCellStore.cpp
	mex FindCentroidDppValue.cpp
	mex Stiff2D.cpp
	mex Stiff2DMonly.cpp
	mex UpdateFEMNodePositions.cpp
	mex UpdatePos.cpp
	mex MeanMaxMinStdTotal.cpp
	mex T1Swaps.cpp
    mex UpdatePosFarhadifar.cpp

end