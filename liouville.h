#ifndef LIOUVILLE_H
#define LIOUVILLE_H

#include <vector>
#include <string>

#include <boost/random/mersenne_twister.hpp>

#include "utilities.h"

#include "Array.h"
#include "fftw++.h"

#define PI 3.141592653589
#define TWOBYPI 0.636619772367581

class RNGenerator {
public:
	RNGenerator();
	double RandomNormal();
	double RandomReal(double min, double max);
	double RandomExponential();
	Complex RandomComplex();
	void Seed(unsigned int seed);
private:
	boost::mt19937 rng_;
};


class GaussianField {
public:
	GaussianField(int width);

	void GenerateRandom(RNGenerator & rng);

	double getField(int x, int y) const;
	int getWidth() const;

	void PrecomputeRandom(RNGenerator & rng);
	void SetMaxMomentum(double momentum);

	void GetFieldList(std::ostream & stream) const;
	void OutputFieldHistogram(std::ostream & stream,int i);
	void DoFieldMeasurement(int i);

private:
	Array::array2<Complex> field_;
	int width_;
	bool max_momentum_changed_;
	bool usereal_;
	void FillSpectrum();
	std::vector<std::vector<double> > spectrum_;
	std::vector<std::vector<double> > measure_;
	std::vector<std::vector<double> > sqrtmeasure_;
	std::vector<double> cumul_;

	std::vector<int> field_measurements_;
	std::vector<std::vector<int> > field_histogram_;
	double min_field_, max_field_;
	int field_bins_;

	bool use_precomputed_random_;
	std::vector<std::vector<Complex> > random_normal_;
	bool use_max_momentum_;
	double max_momentum_;

	fftwpp::fft2d fourier_;
};

class VolumeMeasure {
public:
	VolumeMeasure(int width);
	void setFromGaussianField(const GaussianField & field, double gamma);
	void NormalizeToUnitVolume();

	enum AveragingMethod {
		BOX_SUBDIVISION,
		BOX_AVERAGING,
		DISK_AVERAGING,
		DISK_DIFFUSION
	};

	void performAveraging( VolumeMeasure & target, AveragingMethod method, double delta );
		
	std::pair<int,int> RandomInMeasure(RNGenerator & rng);

	double getMeasure(int x, int y) const;
	double getSqrtMeasure(int x, int y) const;
	void setMeasure(int x, int y, double measure);
	void addMeasure(int x, int y, double measure);
	void updateSqrtMeasure();
		
	int getWidth() const;
	double getTotalVolume() const;
	void Flip();
	void setToZero();
private:
	int width_;

	std::vector<std::vector<double> > measure_;
	std::vector<std::vector<double> > sqrtmeasure_;

	bool cumul_uptodate_;
	std::vector<double> cumul_;
	void UpdateCumulativeMeasure();

	void UpdateSummedMeasure();
	double GetSummedMeasure(int x0, int y0) const;
	double GetSummedMeasure(int x0, int y0, int x1, int y1) const;
	std::vector<std::vector<double> > summed_measure_;
	bool summed_measure_uptodate_;

	std::vector<std::pair<int,int> > spiral_;
	void PrepareSpiral();
	void doBoxSubdivision(VolumeMeasure & target, double delta);
	void doBoxSubdivision(VolumeMeasure & target, int x0, int x1, int y0, int y1, double delta);
	void doBoxAveraging(VolumeMeasure & target, double delta);
	void doDiskAveraging(VolumeMeasure & target, double delta);
	void doDiskDiffusion(VolumeMeasure & target, double delta);
};

struct Node
{
	double distance;
	int x[2];
	short omega[2];
	bool operator<(const Node & node) const { return distance > node.distance; }
};

typedef boost::array<int,2> Vertex;
struct Edge {
	boost::array<Vertex,2> x; 
	bool belongsToLeftFace;
};

class ShortestCycle {
public:
	enum DISTANCE_METHOD {
		DIJKSTRA,
		EIKONAL
	};

	ShortestCycle(VolumeMeasure * const measure, DISTANCE_METHOD method=EIKONAL,double maxLength = 1.5);
	double FindLength();
	void RetrieveCycle(std::vector<Vertex> & path);
	void SetMethod(DISTANCE_METHOD method);
	std::string Output();
private:
	double FindLengthDijkstra();
	double FindLengthEikonal();
	double FindLengthDijkstra(int startX, int startY); 
	double FindLengthEikonal(int startX, int startY); 
	void HorizontalEikonalDistance(int startX);


	std::pair<int,int> Neighbour(const std::pair<int,int> & x, int dir);
	std::pair<short,short> NeighbourOmega( std::pair<int,int> x1, std::pair<int,int> x2, std::pair<short,short> omega ) const;
	bool CrossesBoundary(int boundaryX, int n1x, int n2x) const;

	VolumeMeasure * const measure_;

	std::vector<std::vector<double> > distance_;
	int width_;
	DISTANCE_METHOD method_;

	std::vector<int> histogram_;
	double average_length_;
	double min_length_, max_length_;
	int bins_;
	int measurements_;

	Vertex cycle_base_point_;
	boost::array<Vertex,2> cycle_tips_;
	boost::array<Edge,2 > cycle_edges_;
};

class EikonalGeodesicRetriever {
public:
	EikonalGeodesicRetriever(std::vector<std::vector<double> > & distance);
	void FindGeodesicFromTips(std::vector<Vertex> & path, const boost::array<Vertex,2> & tips);
	void FindGeodesic(std::vector<Vertex> & path);
	void FindGeodesic(std::vector<Vertex> & path, Edge edge);
private:
	std::vector<std::vector<double> > & distance_;
	int width_;

	void nextEdge(Edge & edge,double & lambda) const;
	Edge Adjacent(const Edge & edge) const;
	Edge Next(const Edge & edge) const;
	Edge Previous(const Edge & edge) const;
	Edge Reverse(const Edge & edge) const;
	double DistanceChange(const Edge & edge) const;
};

class TwoPointFunction {
public:
	enum DISTANCE_METHOD {
		DIJKSTRA,
		EIKONAL
	};

	TwoPointFunction(VolumeMeasure * const measure, DISTANCE_METHOD method, double maxdist = 1.3, int bins = 1000);
	void Measure(RNGenerator & rng, int n=1);

	double getDistance(int x, int y) const;
	std::string Output();
	void DijkstraDistance(int startX, int startY, bool measurement=true);
	void EikonalDistance(int startX, int startY, bool measurement=true);

	void SetMethod(DISTANCE_METHOD method);
	void SaveDistance(std::string filename) const;
private:
	VolumeMeasure * const measure_;
	std::pair<int,int> Neighbour(const std::pair<int,int> & x, int dir);

	std::vector<std::vector<double> > distance_;

	int width_;

	DISTANCE_METHOD method_;

	double max_dist_;
	int bins_;
	std::vector<double> dist_histogram_;
	int measurements_;

};

class CorrelationFunction {
public:
	CorrelationFunction(VolumeMeasure * const measure);
	void MeasureCorrelation(RNGenerator & rng, int n=1);
	std::string Output();

private:
	std::vector<double> correlation_;
	int correlation_measurements_;
	VolumeMeasure * const measure_;
	int width_;
};


#endif