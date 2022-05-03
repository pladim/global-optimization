#ifndef HELPFUL_STRUCTS_H
#define HELPFUL_STRUCTS_H

struct SummaryPacket {
	wchar_t* summary;
	int length;
};

struct GKLSGParametersPacket {
	int numberProblem;
	int dimention;
	int maxNumberConstraints;
	int numberLocalMinimasMinimandFunction;
	int numberLocalMinimasConstraintFunction;
	int minNumberActiveConstraints;
	int maxNumberActiveConstraints;
	double distanceFromGlobalMinimizerToParaboloidMinimizerMinimandFunction;
	double minDistanceFromMinimizerToParaboloidMinimizer;
	double maxDistanceFromMinimizerToParaboloidMinimizer;
	double radiusAttractionRegionGlobalMinimizerMminimandFunction;
	double minRadiusAttractionRegionGlobalMinimizer;
	double maxRadiusAttractionRegionGlobalMinimizer;
	double valueGlobalMinimumMinimandFunction;
	double minValueGlobalMinimumConstraints;
	double maxValueGlobalMinimumConstraints;
	double minAnglePhi;
	double maxAnglePhi;
	double maxPower;
	double probabilityConvexityConstraints;
	double minNu0;
	double maxNu0;
	double minDeltaGlobal;
	double maxDeltaGlobal;
	double deltaLocal;
	double maxLambda;
	double probabilityGglobalMinimumIsOutsidePermissibleRegion;
	double relativeMeasurePermissibleRegion;
	int functionType;
};

#endif // HELPFUL_STRUCTS_H
