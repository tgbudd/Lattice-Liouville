#include <time.h>
#include <algorithm>
#include <vector>
#include <queue>
#include <string>
#include <sstream>
#include <fstream>

#include <boost/assert.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/uniform_real.hpp>

#include "liouville.h"
#include "utilities.h"

RNGenerator::RNGenerator() :
	rng_(static_cast<unsigned int>(time(NULL)))
{
}

double RNGenerator::RandomNormal()
{
  	boost::normal_distribution<> normal_distribution(0.0,1.0);
	boost::variate_generator< boost::mt19937&, boost::normal_distribution<> > rnGenerator(rng_,normal_distribution);
	double random = rnGenerator();
	return random;
}

Complex RNGenerator::RandomComplex()
{
	return Complex(RandomNormal(),RandomNormal());
}

double RNGenerator::RandomReal(double min, double max)
{
	boost::uniform_real<> distribution(min,max);
	return distribution(rng_);
}

double RNGenerator::RandomExponential()
{
	return -std::log(RandomReal(0.0000001,1.0));
}

void RNGenerator::Seed(unsigned int seed)
{
	rng_.seed(seed);
}


GaussianField::GaussianField(int width) :
	field_(width,width,sizeof(Complex)),
	width_(width),
	spectrum_(width_,std::vector<double>(width_,0.0)),
	usereal_(false),
	fourier_(1,field_),
	use_precomputed_random_(false),
	use_max_momentum_(false),
	max_momentum_changed_(true),
	min_field_(-10.0),
	max_field_(10.0),
	field_bins_(800)
{
	FillSpectrum();
}

void GaussianField::GenerateRandom(RNGenerator & rng)
{
	if( usereal_ && !(use_max_momentum_ && max_momentum_changed_) )
	{
		// Use the imaginary part of the field.
		usereal_ = false;
	} else
	{
		// Generate a new complex gaussian field
		// and use its real part.
		usereal_ = true;
		double rescaling = 1.0/width_;
		for(int p1=0;p1<width_;p1++)
		{
			for(int p2=0;p2<width_;p2++)
			{
				if( use_max_momentum_ && max_momentum_ * spectrum_[p1][p2] < 1.0 )
				{
					field_(p1,p2) = 0.0;
				} else
				{
					field_(p1,p2) = rescaling * spectrum_[p1][p2] 
						* (use_precomputed_random_ ? random_normal_[p1][p2] : rng.RandomComplex() );
				}
			}
		}
		fourier_.fft(field_);
		max_momentum_changed_ = false;
	}
}

double GaussianField::getField(int x, int y) const
{
	return (usereal_ ? field_(x,y).real() : field_(x,y).imag());
}

int GaussianField::getWidth() const
{
	return width_;
}

void GaussianField::SetMaxMomentum(double momentum)
{
	use_max_momentum_ = true;
	max_momentum_ = momentum;
	max_momentum_changed_ = true;
}

void GaussianField::FillSpectrum()
{
	std::vector<double> sines(width_);
	for(int p=0;p<width_;p++)
	{
		double sine = std::sin(static_cast<double>(p) * PI / static_cast<double>(width_));
		sines[p] = TWOBYPI * sine * sine;
	}

	for(int p1=0;p1<width_;p1++)
	{
		for(int p2=0;p2<width_;p2++)
		{
			if( p1 != 0 || p2 != 0 )
				spectrum_[p1][p2] = 1.0/std::sqrt( sines[p1] + sines[p2] );
			else
				spectrum_[p1][p2] = 0.0;
		}
	}
}


void GaussianField::PrecomputeRandom(RNGenerator & rng)
{
	use_precomputed_random_ = true;
	if(	static_cast<int>(random_normal_.size()) != width_ )
	{
		random_normal_.resize(width_,std::vector<Complex>(width_));
	}
	for(int i=0;i<width_;i++)
	{
		for(int j=0;j<width_;j++)
		{
			random_normal_[i][j] = rng.RandomComplex();
		}
	}
}

void GaussianField::GetFieldList(std::ostream & stream) const
{
	stream << "{" << std::fixed;
	for(int i=0;i<width_;i++)
	{
		stream << (i>0?",":"") << "{";
		for(int j=0;j<width_;j++)
		{
			stream << (j>0?",":"") << (usereal_?field_(i,j).real() : field_(i,j).imag());
		}
		stream << "}";
	}
	stream << "}";
}


void GaussianField::DoFieldMeasurement(int index)
{
	if( index >= static_cast<int>(field_measurements_.size()) )
	{
		field_measurements_.resize(index+1,0);
		field_histogram_.resize(index+1,std::vector<int>(field_bins_,0));
	}

	for(int x=0;x<width_;x++)
	{
		for(int y=0;y<width_;y++)
		{
			double field = getField(x,y);
			int n = static_cast<int>((field-min_field_)/(max_field_-min_field_)*field_bins_);
			if( n >= 0 && n < field_bins_ )
			{
				field_histogram_[index][n]++;
			}
		}
	}

	field_measurements_[index]++;
}

void GaussianField::OutputFieldHistogram(std::ostream & stream,int i)
{
	stream << std::fixed << "{ measurements -> " << field_measurements_[i] << ", minfield -> " << min_field_;
	stream << ", maxfield -> " << max_field_ << ", fieldbins -> " << field_bins_ << ", fieldhistogram -> ";
	PrintToStream(stream,field_histogram_[i].begin(),field_histogram_[i].end());
	stream << "}";
}


VolumeMeasure::VolumeMeasure(int width) :
	width_(width),
	summed_measure_uptodate_(false),
	cumul_uptodate_(false)
{
	measure_.resize(width_,std::vector<double>(width_));
	sqrtmeasure_.resize(width_,std::vector<double>(width_));
	cumul_.resize(width_*width_);
}

double VolumeMeasure::getMeasure(int x, int y) const
{
	return measure_[x][y];
}
double VolumeMeasure::getSqrtMeasure(int x, int y) const
{
	return sqrtmeasure_[x][y];
}

int VolumeMeasure::getWidth() const
{
	return width_;
}

void VolumeMeasure::setFromGaussianField(const GaussianField & field, double gamma)
{
	// Set the measure according to a gaussian free field.
	// Notice that measure_[x][y] corresponds to the Liouville measure integrated
	// over the lattice site at position (x,y), i.e. it is NOT a density.
	double renorm = std::pow(static_cast<double>(width_),-2.0-gamma*gamma/2.0);  

	for(int p1=0;p1<width_;p1++)
	{
		for(int p2=0;p2<width_;p2++)
		{
			double measure = renorm * std::exp(gamma * field.getField(p1,p2) );
			setMeasure(p1,p2, measure );
		}
	}
}

void VolumeMeasure::setMeasure(int x, int y, double measure)
{
	measure_[x][y] = measure;
	sqrtmeasure_[x][y] = std::sqrt(measure);
	summed_measure_uptodate_ = false;
	cumul_uptodate_ = false;
}

void VolumeMeasure::addMeasure(int x, int y, double measure)
{
	measure_[x][y] += measure;
	summed_measure_uptodate_ = false;
	cumul_uptodate_ = false;
}

void VolumeMeasure::NormalizeToUnitVolume()
{
	double volume = getTotalVolume();
	for(int x=0;x<width_;x++)
	{
		for(int y=0;y<width_;y++)
		{
			setMeasure(x,y,getMeasure(x,y)/volume);
		}
	}
}

void VolumeMeasure::updateSqrtMeasure()
{
	for(int x=0;x<width_;x++)
	{
		for(int y=0;y<width_;y++)
		{
			sqrtmeasure_[x][y] = std::sqrt(measure_[x][y]);
		}
	}
}

void VolumeMeasure::UpdateCumulativeMeasure()
{
	double total = 0.0;
	for(int p1=0;p1<width_;p1++)
	{
		for(int p2=0;p2<width_;p2++)
		{
			total += getMeasure(p1,p2);
			cumul_[p2 + p1 * width_] = total;
		}
	}
	cumul_uptodate_ = true;
}

std::pair<int,int> VolumeMeasure::RandomInMeasure(RNGenerator & rng)
{
	if( !cumul_uptodate_ )
	{
		UpdateCumulativeMeasure();
	}
	double x = rng.RandomReal(0.0,cumul_.back());

	// Perform a binary search to locate the first lattice site 
	// for which the cumulative measure exceeds x.
	int imin = 0, imax = cumul_.size()-1;
	double min = cumul_[imin], max = cumul_[imax];
	while( min < x )
	{
		int next = static_cast<int>(imin + (imax-imin)*((x-min)/(max-min)));

		if( cumul_[next] < x )
		{
			imin = next + 1;
			min = cumul_[imin];
		} else
		{
			imax = next;
			max = cumul_[imax];
		}
	}
	BOOST_ASSERT( cumul_[imin] >= x );
	BOOST_ASSERT( imin == 0 || cumul_[imin-1] < x );
	return std::pair<int,int>(imin / width_, imin % width_);
}

void VolumeMeasure::performAveraging( VolumeMeasure & target, AveragingMethod method, double delta )
{
	if( method == BOX_SUBDIVISION )
	{
		doBoxSubdivision(target,delta);
	} else if( method == BOX_AVERAGING )
	{
		doBoxAveraging(target,delta);
	} else if( method == DISK_AVERAGING )
	{
		doDiskAveraging(target,delta);
	} else if( method == DISK_DIFFUSION )
	{
		doDiskDiffusion(target,delta);
	}
}

void VolumeMeasure::doBoxSubdivision( VolumeMeasure & target, double delta)
{
	bool powerOfTwo = !(width_ == 0) && !(width_ & (width_ - 1));
	BOOST_ASSERT( powerOfTwo );

	if( !summed_measure_uptodate_ )
	{
		UpdateSummedMeasure();
	}

	doBoxSubdivision(target, 0,width_,0,width_,delta);
}

void VolumeMeasure::doBoxSubdivision( VolumeMeasure & target, int x0, int x1, int y0, int y1, double delta)
{
	if( x1-x0 == 1 )
	{
		target.setMeasure(x0,y0,getMeasure(x0,y0));
	} else
	{
		double measure = GetSummedMeasure(x0,y0,x1-1,y1-1);
		if( measure < delta )
		{
			double averagemeasure = measure/((x1-x0)*(y1-y0));
			BOOST_ASSERT(averagemeasure > 0.0);
			for(int x=x0;x<x1;x++)
			{
				for(int y=y0;y<y1;y++)
				{
					target.setMeasure(x,y,averagemeasure);
				}
			}
			int logsize = logtwo(x1-x0);
		} else
		{
			doBoxSubdivision(target,x0,(x0+x1)/2,y0,(y0+y1)/2,delta);
			doBoxSubdivision(target,x0,(x0+x1)/2,(y0+y1)/2,y1,delta);
			doBoxSubdivision(target,(x0+x1)/2,x1,y0,(y0+y1)/2,delta);
			doBoxSubdivision(target,(x0+x1)/2,x1,(y0+y1)/2,y1,delta);
		}
	}
}

void VolumeMeasure::doBoxAveraging(VolumeMeasure & target, double delta)
{
	if( !summed_measure_uptodate_ )
	{
		UpdateSummedMeasure();
	}

	for(int x=0;x<width_;x++)
	{
		for(int y=0;y<width_;y++)
		{
			double lastmeasure = 0.0;
			for(int r=0;r<width_/2;r++)
			{
				double measure = GetSummedMeasure( x-r,y-r,x+r,y+r );
				if( measure > delta )
				{
					double radius = std::max(0.5,static_cast<double>(r)-0.5+(delta-lastmeasure)/(measure-lastmeasure));
					target.setMeasure(x,y,delta/(4.0*radius*radius)); 
					break;
				}
				lastmeasure = measure;
			}
		}
	}
}

void VolumeMeasure::UpdateSummedMeasure()
{
	if( static_cast<int>(summed_measure_.size()) != width_ )
	{
		summed_measure_.resize(width_,std::vector<double>(width_,0.0));
	}
	for(int x=0;x<width_;x++)
	{
		for(int y=0;y<width_;y++)
		{
			summed_measure_[x][y] = getMeasure(x,y) + GetSummedMeasure(x-1,y) + GetSummedMeasure(x,y-1) - GetSummedMeasure(x-1,y-1);
		}
	}
	summed_measure_uptodate_ = true;
}

double VolumeMeasure::GetSummedMeasure(int x0,int y0) const
{
	return (x0 < 0 || y0 < 0 || x0 >= width_ || y0 >= width_ ? 0.0 : summed_measure_[x0][y0] );
}

double VolumeMeasure::GetSummedMeasure(int x0, int y0, int x1, int y1) const
{
	if( x0 < 0 || y0 < 0 || x0 >= width_ || y0 >= width_ )
	{
		int shiftx = ProperMod(x0,width_)-x0;
		int shifty = ProperMod(y0,width_)-y0;
		return GetSummedMeasure( x0 + shiftx, y0 + shifty, x1 + shiftx, y1 + shifty );
	}

	if( x1 >= width_ )
	{
		return GetSummedMeasure( x0, y0, width_-1, y1 ) + GetSummedMeasure( 0, y0, x1-width_, y1 );
	}

	if( y1 >= width_ )
	{
		return GetSummedMeasure( x0, y0, x1, width_-1 ) + GetSummedMeasure( x0, 0, x1, y1-width_ );
	}

	return GetSummedMeasure(x1,y1) + GetSummedMeasure(x0-1, y0-1) - GetSummedMeasure(x1,y0-1) - GetSummedMeasure(x0-1,y1);
}



void VolumeMeasure::doDiskAveraging(VolumeMeasure & target, double delta)
{
	if( spiral_.empty() )
	{
		PrepareSpiral();
	}
	for(int x=0;x<width_;x++)
	{
		for(int y=0;y<width_;y++)
		{
			double measure = 0.0;
			for(int a=0,enda=spiral_.size();a<enda;a++)
			{
				measure += getMeasure(ProperMod(x+spiral_[a].first,width_),ProperMod(y+spiral_[a].second,width_));
				if( measure > delta )
				{
					if( a == 0 )
					{
						target.setMeasure(x,y,measure);
					} else
					{
						target.setMeasure(x,y,delta/(a+1));
					}
					//target.setMeasure(x,y,measure/(a+1));
					break;
				}
			}
		}
	}
}

void VolumeMeasure::setToZero()
{
	for(int x=0;x<width_;x++)
	{
		for(int y=0;y<width_;y++)
		{
			setMeasure(x,y,0.0);
		}
	}
}

void VolumeMeasure::doDiskDiffusion(VolumeMeasure & target, double delta)
{
	if( spiral_.empty() )
	{
		PrepareSpiral();
	}
	target.setToZero();
	for(int x=0;x<width_;x++)
	{
		for(int y=0;y<width_;y++)
		{
			double measure = 0.0;
			int ballsize=1;
			for(int a=0,enda=spiral_.size();a<enda;a++)
			{
				measure += getMeasure(ProperMod(x+spiral_[a].first,width_),ProperMod(y+spiral_[a].second,width_));
				if( measure >= delta )
				{
					ballsize = a+1;
					break;
				}
			}
			double add = getMeasure(x,y)/ballsize;
			for(int a=0;a<ballsize;a++)
			{
				target.addMeasure(ProperMod(x+spiral_[a].first,width_),ProperMod(y+spiral_[a].second,width_),add);
			}
		}
	}
	target.updateSqrtMeasure();
}

bool wayToSort(const std::pair<int,std::pair<int,int> > & a, const std::pair<int,std::pair<int,int> > & b)
{
	return a.first < b.first;
}

void VolumeMeasure::PrepareSpiral()
{
	std::vector<std::pair<int,std::pair<int,int> > > dists;
	for(int x=-width_/2;x<width_/2;x++)
	{
		for(int y=-width_/2;y<width_/2;y++)
		{
			int sqdist = x*x + y*y;
			if( sqdist < (width_/2)*(width_/2) )
			{
				dists.push_back(std::pair<int,std::pair<int,int> >(sqdist,std::pair<int,int>(x,y)));
			}
		}
	}
	std::sort(dists.begin(),dists.end(),wayToSort);
	for(int i=0,endi=dists.size();i<endi;i++)
	{
		spiral_.push_back(dists[i].second);
	}
}

void VolumeMeasure::Flip()
{
	// mirror in the diagonal
	for(int x=0;x<width_;x++)
	{
		for(int y=x+1;y<width_;y++)
		{
			std::swap(measure_[x][y],measure_[y][x]);
			std::swap(sqrtmeasure_[x][y],sqrtmeasure_[y][x]);
		}
	}
	summed_measure_uptodate_ = false;
	cumul_uptodate_ = false;
}

double VolumeMeasure::getTotalVolume() const
{
	double total=0.0;
	for(int x=0;x<width_;x++)
	{
		for(int y=0;y<width_;y++)
		{
			total+=measure_[x][y];
		}
	}
	return total;
}

ShortestCycle::ShortestCycle(VolumeMeasure * const measure, DISTANCE_METHOD method, double maxLength) :
	measure_(measure),
	method_(method),
	distance_(measure->getWidth(),std::vector<double>(measure->getWidth())),
	width_(measure->getWidth()),
	measurements_(0),
	average_length_(0.0),
	min_length_(0.0),
	max_length_(maxLength),
	bins_(1000)
{
	histogram_.resize(bins_,0);
}

void ShortestCycle::SetMethod(DISTANCE_METHOD method)
{
	method_ = method;
}

double ShortestCycle::FindLength()
{
	double shortestcycle;
	if( method_ == DIJKSTRA )
	{
		shortestcycle = FindLengthDijkstra();
	} else
	{
		shortestcycle = FindLengthEikonal();
	}

	average_length_ += shortestcycle;
	int bin = static_cast<int>(bins_*(shortestcycle - min_length_)/(max_length_ - min_length_));
	if( bin >=0 && bin < bins_ )
	{
		histogram_[bin]++;
	}
	measurements_++;

	return shortestcycle;
}

double ShortestCycle::FindLengthDijkstra(int startX, int startY)
{
	const double VERYLARGE = 10000.0;
	for(int i=0;i<width_;i++)
	{
		for(int j=0;j<width_;j++)
		{
			distance_[i][j] = VERYLARGE;
		}
	}

	std::priority_queue<Node> pqueue;
	Node startNode;
	startNode.distance = 0.0;
	startNode.x[0] = startX;
	startNode.x[1] = startY;
	startNode.omega[0] = 0;
	startNode.omega[1] = 0;
	pqueue.push( startNode );
	distance_[startX][startY] = 0.0;
	std::vector<std::vector<std::pair<short,short> > > omega(width_,std::vector<std::pair<short,short> >(width_));
	omega[startX][startY].first = 0;
	omega[startX][startY].second = 0;
	double shortestcycle = VERYLARGE;
	boost::array<boost::array<int,2>,2 > cyclenodes;

	while( !pqueue.empty() )
	{
		Node node = pqueue.top();
		pqueue.pop();

		if( node.distance > distance_[node.x[0]][node.x[1]] )
			continue;

		for(int i=0;i<4;i++)
		{
			Node nextNode(node);
			if( i==0 )
				nextNode.x[0] = (nextNode.x[0]+1)%width_;
			else if( i==1 )
				nextNode.x[0] = (nextNode.x[0]+width_-1)%width_;
			else if( i==2 )
				nextNode.x[1] = (nextNode.x[1]+1)%width_;
			else
				nextNode.x[1] = (nextNode.x[1]+width_-1)%width_;
			for(int j=0;j<2;j++)
			{
				if( node.x[j] == 0 && nextNode.x[j] == width_-1 )
					nextNode.omega[j]--;
				else if( node.x[j] == width_-1 && nextNode.x[j] == 0 )
					nextNode.omega[j]++;
			}

			nextNode.distance += 0.5 * ( measure_->getSqrtMeasure(node.x[0],node.x[1]) + measure_->getSqrtMeasure(nextNode.x[0],nextNode.x[1]) );

			if( nextNode.distance < distance_[nextNode.x[0]][nextNode.x[1]] )
			{
				distance_[nextNode.x[0]][nextNode.x[1]] = nextNode.distance;
				omega[nextNode.x[0]][nextNode.x[1]].first = nextNode.omega[0];
				omega[nextNode.x[0]][nextNode.x[1]].second = nextNode.omega[1];
				pqueue.push( nextNode );
			} else if ( distance_[nextNode.x[0]][nextNode.x[1]] <= node.distance      // the distance is already fixed ...
				&& (nextNode.omega[0] != omega[nextNode.x[0]][nextNode.x[1]].first    // ... and originates from a non-trivial cycle
				|| nextNode.omega[1] != omega[nextNode.x[0]][nextNode.x[1]].second ) )
			{
				double length = nextNode.distance + distance_[nextNode.x[0]][nextNode.x[1]];
				if( length < shortestcycle )
				{
					shortestcycle = length;
					cyclenodes[0][0] = node.x[0];
					cyclenodes[0][1] = node.x[1];
					cyclenodes[1][0] = nextNode.x[0];
					cyclenodes[1][1] = nextNode.x[1];
				}
			}
		}
	}
	return shortestcycle;
}

double ShortestCycle::FindLengthDijkstra()
{
	double shortestcycle = 10000.0;
	std::vector<double> tmp;
	for(int y=0;y<width_;y++)
	{
		double cycle = FindLengthDijkstra(width_-1,y);
		tmp.push_back(cycle);
		if( cycle < shortestcycle )
		{
			shortestcycle = cycle;
		}
	}
	for(int x=0;x<width_-1;x++)
	{
		double cycle = FindLengthDijkstra(x,width_-1);
		tmp.push_back(cycle);
		if( cycle < shortestcycle )
		{
			shortestcycle = cycle;
		}
	}
	return shortestcycle;
}

template <class T>
class OrderByVector {
public:
	OrderByVector(const std::vector<T> & v) : v_(v) {}
	bool operator() (int i, int j) {
		return v_[i] < v_[j];
	}
private:
	std::vector<T> v_;
};

double ShortestCycle::FindLengthEikonal()
{
	double shortestcycle = 10000.0;
	std::vector<double> lowerbound;
	std::vector<Vertex> points(2*width_-1);
	std::vector<int> start(2*width_-1);
	for(int y=0;y<2*width_-1;y++)
	{
		start[y]=y;
		points[y][0]=(y<width_?width_-1:y-width_+1);
		points[y][1]=(y<width_?y:width_-1);
	}
	for(int j=0;j<2;j++)
	{
		HorizontalEikonalDistance(0);
		for(int y=j;y<width_;y++)
		{
			lowerbound.push_back(distance_[width_-1][y]);
		}
		measure_->Flip();
	}

	OrderByVector<double> order(lowerbound);
	std::sort(start.begin(),start.end(),order);
	for(int i=0;i<width_;i+=4)
	{
		int y = start[i];
		if( lowerbound[y] < shortestcycle )
		{
			double cycle = FindLengthEikonal(points[y][0],points[y][1]);
			if( cycle < shortestcycle )
			{
				shortestcycle = cycle;
				cycle_base_point_ = points[y];
			}
		}
	}

	return shortestcycle;
}


void ShortestCycle::HorizontalEikonalDistance(int startX)
{
	for(int i=0;i<width_;i++)
	{
		for(int j=0;j<width_;j++)
		{
			distance_[i][j] = 100000.0;
		}
	}

	std::priority_queue<Node> pqueue;
	std::vector<Node> startNodes(width_);
	for(int i=0;i<width_;i++)
	{
		startNodes[i].x[0] = startX;
		startNodes[i].x[1] = i;
		startNodes[i].distance = 0.0;
		pqueue.push( startNodes[i] );
		distance_[startX][i] = 0.0;
	}
	std::vector<double> mindistances(width_,10000.0);

	while( !pqueue.empty() )
	{
		Node node = pqueue.top();
		pqueue.pop();

		if( node.distance > distance_[node.x[0]][node.x[1]] )
			continue;

		if( node.distance < mindistances[node.x[0]] )
		{
			mindistances[node.x[0]] = node.distance;
		}

		for(int i=0;i<4;i++)
		{
			Node nextNode(node);
			std::pair<int,int> nextx = Neighbour(std::pair<int,int>(node.x[0],node.x[1]),i);
			nextNode.x[0] = nextx.first;
			nextNode.x[1] = nextx.second;

			if( CrossesBoundary(startX,node.x[0],nextNode.x[0] ) )
			{
				continue;
			}

			if( distance_[nextNode.x[0]][nextNode.x[1]] > node.distance )
			{
				// its distance is not fixed yet
				std::pair<int,int> nextx2 = Neighbour(nextx,i);
				double a = std::min( node.distance, ( CrossesBoundary(startX,nextNode.x[0],nextx2.first) ? 100000.0 : distance_[nextx2.first][nextx2.second] ) );
				std::pair<int,int> nextx3 = Neighbour(nextx,2*(1-i/2));
				std::pair<int,int> nextx4 = Neighbour(nextx,2*(1-i/2)+1);
				double b = std::min( ( CrossesBoundary(startX,nextNode.x[0],nextx3.first) ? 100000.0 : distance_[nextx3.first][nextx3.second] ), 
					( CrossesBoundary(startX,nextNode.x[0],nextx4.first) ? 100000.0 : distance_[nextx4.first][nextx4.second] ) );
				
				double localDistance = measure_->getSqrtMeasure(nextNode.x[0],nextNode.x[1]);
				if( a > 10000.0 || b > 10000.0 )
				{
					nextNode.distance = std::min(a,b) + localDistance;
				} else
				{
					double discriminant = 2 * localDistance * localDistance - (a-b)*(a-b);
					if( discriminant >= 0.0 )
					{
						nextNode.distance = 0.5*(a+b+std::sqrt(discriminant));
					} else
					{
						nextNode.distance = std::min(a,b) + localDistance;
					}
				}


				if( nextNode.distance < distance_[nextNode.x[0]][nextNode.x[1]] )
				{
					distance_[nextNode.x[0]][nextNode.x[1]] = nextNode.distance;
					pqueue.push( nextNode );
				}
			}
		}
	}
}

void ShortestCycle::RetrieveCycle(std::vector<Vertex> & path)
{
	// For now only eikonal geodesics can be retrieved
	BOOST_ASSERT( method_ == EIKONAL );

	FindLengthEikonal(cycle_base_point_[0],cycle_base_point_[1]);

	EikonalGeodesicRetriever geodesic(distance_);
	geodesic.FindGeodesicFromTips(path,cycle_tips_);
}

std::pair<short,short> ShortestCycle::NeighbourOmega( std::pair<int,int> x1, std::pair<int,int> x2, std::pair<short,short> omega ) const
{
	if( x1.first == 0 && x2.first == width_-1 )
	{
		omega.first--;
	} else if( x1.first == width_-1 && x2.first == 0 )
	{
		omega.first++;
	}
	if( x1.second == 0 && x2.second == width_-1 )
	{
		omega.second--;
	} else if( x1.second == width_-1 && x2.second == 0 )
	{
		omega.second++;
	}
	return omega;
}

double ShortestCycle::FindLengthEikonal(int startX, int startY)
{
	const double VERYLARGE = 10000.0;
	for(int i=0;i<width_;i++)
	{
		for(int j=0;j<width_;j++)
		{
			distance_[i][j] = VERYLARGE+1.0;
		}
	}
	std::priority_queue<Node> pqueue;
	Node startNode;
	startNode.distance = 0.0;
	startNode.x[0] = startX;
	startNode.x[1] = startY;
	startNode.omega[0] = 0;
	startNode.omega[1] = 0;
	pqueue.push( startNode );
	distance_[startX][startY] = 0.0;
	std::vector<std::vector<std::pair<short,short> > > omega(width_,std::vector<std::pair<short,short> >(width_));
	omega[startX][startY].first = 0;
	omega[startX][startY].second = 0;
	double shortestcycle = VERYLARGE;
	double current_distance_ = 0.0;

	while( !pqueue.empty() )
	{
		Node node = pqueue.top();
		pqueue.pop();

		if( node.distance > distance_[node.x[0]][node.x[1]] )
			continue;

		current_distance_ = node.distance;

		for(int i=0;i<4;i++)
		{
			Node nextNode(node);
			std::pair<int,int> nextx = Neighbour(std::pair<int,int>(node.x[0],node.x[1]),i);
			nextNode.x[0] = nextx.first;
			nextNode.x[1] = nextx.second;
			std::pair<short,short> nextomega = NeighbourOmega(std::pair<int,int>(node.x[0],node.x[1]),nextx,std::pair<short,short>(node.omega[0],node.omega[1]));
			nextNode.omega[0] = nextomega.first;
			nextNode.omega[1] = nextomega.second;
	
			if( distance_[nextNode.x[0]][nextNode.x[1]] > current_distance_ || nextNode.omega[0] != omega[nextNode.x[0]][nextNode.x[1]].first 
					|| nextNode.omega[1] != omega[nextNode.x[0]][nextNode.x[1]].second  )
			{
				// its distance is not fixed yet
				std::pair<int,int> nextx2 = Neighbour(nextx,i);
				std::pair<short,short> nextomega2 = NeighbourOmega(nextx,nextx2,nextomega);
				double a = std::min( node.distance, (nextomega2 == omega[nextx2.first][nextx2.second] ? distance_[nextx2.first][nextx2.second] : VERYLARGE ) );
				std::pair<int,int> nextx3 = Neighbour(nextx,2*(1-i/2));
				std::pair<short,short> nextomega3 = NeighbourOmega(nextx,nextx3,nextomega);
				std::pair<int,int> nextx4 = Neighbour(nextx,2*(1-i/2)+1);
				std::pair<short,short> nextomega4 = NeighbourOmega(nextx,nextx4,nextomega);
				double b = std::min( (nextomega3 == omega[nextx3.first][nextx3.second] ? distance_[nextx3.first][nextx3.second] : VERYLARGE ),
					                 (nextomega4 == omega[nextx4.first][nextx4.second] ? distance_[nextx4.first][nextx4.second] : VERYLARGE ) );
				
				double localDistance = measure_->getSqrtMeasure(nextNode.x[0],nextNode.x[1]);
				if( a > 10000.0 || b > 10000.0 )
				{
					nextNode.distance = std::min(a,b) + localDistance;
				} else
				{
					double discriminant = 2 * localDistance * localDistance - (a-b)*(a-b);
					if( discriminant >= 0.0 )
					{
						nextNode.distance = 0.5*(a+b+std::sqrt(discriminant));
					} else
					{
						nextNode.distance = std::min(a,b) + localDistance;
					}
				}

				if( nextNode.distance < distance_[nextNode.x[0]][nextNode.x[1]] )
				{
					distance_[nextNode.x[0]][nextNode.x[1]] = nextNode.distance;
					omega[nextNode.x[0]][nextNode.x[1]].first = nextNode.omega[0];
					omega[nextNode.x[0]][nextNode.x[1]].second = nextNode.omega[1];
					pqueue.push( nextNode );
				} else if ( distance_[nextNode.x[0]][nextNode.x[1]] <= current_distance_ )
				{
					double length = nextNode.distance + distance_[nextNode.x[0]][nextNode.x[1]];
					if( length < shortestcycle )
					{
						shortestcycle = length;
						cycle_tips_[0][0] = node.x[0];
						cycle_tips_[0][1] = node.x[1];
						cycle_tips_[1][0] = nextNode.x[0];
						cycle_tips_[1][1] = nextNode.x[1];

					}
				}
			}
		}
	}
	return shortestcycle;
}

bool ShortestCycle::CrossesBoundary(int boundaryX, int n1x, int n2x) const
{
	return (n1x == boundaryX && n2x == ProperMod(boundaryX-1,width_)) || (n2x == boundaryX && n1x == ProperMod(boundaryX-1,width_));
}

std::pair<int,int> ShortestCycle::Neighbour(const std::pair<int,int> & x, int dir)
{
	return std::pair<int,int>( (dir > 1 ? x.first : (dir == 0 ? (x.first+1)%width_ : (x.first+width_-1)%width_ ) ),
		(dir <= 1 ? x.second : (dir == 2 ? (x.second+1)%width_ : (x.second+width_-1)%width_ ) ) );
}

std::string ShortestCycle::Output()
{
	std::ostringstream os;
	os << std::fixed << "{width -> " << width_;
	os << ", method -> \"" << (method_==EIKONAL ? "Eikonal" : "Dijkstra" ) << "\"";
	os << ", measurements -> " << measurements_;
	os << ", minlength -> " << min_length_;
	os << ", maxlength -> " << max_length_;
	os << ", bins -> " << bins_;
	os << ", cyclehistogram -> ";
	PrintToStream(os,histogram_.begin(),histogram_.end());
	os << ", averagecyclelength -> " << average_length_ / measurements_;
	os << ", cyclemeasurements -> " << measurements_;
	os << "}";
	return os.str();
}


EikonalGeodesicRetriever::EikonalGeodesicRetriever(std::vector<std::vector<double> > & distance) :
	distance_(distance),
	width_(distance.size())
{
	
}
	
void EikonalGeodesicRetriever::FindGeodesicFromTips(std::vector<Vertex> & path, const boost::array<Vertex,2> & tips)
{
	Edge edge;
	edge.x = tips;
	edge.belongsToLeftFace = true;

	Edge edge1 = Reverse(Adjacent(Previous(edge)));
	Edge edge2 = Adjacent(Next(Adjacent(edge)));

	path.clear();
	path.push_back(edge1.x[0]);
	FindGeodesic(path,edge1);

	std::vector<Vertex> path2;
	path2.push_back(edge2.x[0]);
	FindGeodesic(path2,edge2);
	std::reverse(path2.begin(),path2.end());
	path.insert(path.begin(),path2.begin(),path2.end());
}

void EikonalGeodesicRetriever::FindGeodesic(std::vector<Vertex> & path)
{
	Edge edge;
	Vertex p = path.back();
	edge.x[0] = p;
	edge.x[1][0] = ProperMod(p[0]+1,width_);
	edge.x[1][1] = p[1];
	edge.belongsToLeftFace = true;
	FindGeodesic(path,edge);
}

void EikonalGeodesicRetriever::FindGeodesic(std::vector<Vertex> & path, Edge edge)
{
	BOOST_ASSERT(!path.empty());
	Vertex p = path.back();
	Vertex lastp = p;
	double lambda = 0.45;
	while( !( distance_[edge.x[0][0]][edge.x[0][1]] < 0.00001 ) && !( distance_[edge.x[1][0]][edge.x[1][1]] < 0.00001 ) )
	{
		nextEdge(edge,lambda);

		if( lambda < 0.5 )
		{
			p = edge.x[0];
		} else
		{
			p = edge.x[1];
		}
		path.push_back(p);
		lastp = p;
	}
}

void EikonalGeodesicRetriever::nextEdge(Edge & edge, double & lambda) const
{
	// Make sure the edge is oriented towards the diagonal.
	// The diagonal is defined to avoid the vertex of minimal distance.
	double mindistance = 0.0, currentdistance = 0.0;
	int minpos = 0;
	for(int i=1;i<=4;i++)
	{
		currentdistance += DistanceChange(edge);
		if( currentdistance < mindistance )
		{
			mindistance = currentdistance;
			minpos = i;
		}
		edge = Next(edge);
	}
	if( minpos == 1 || minpos == 3 )
	{
		edge = Reverse(edge);
		lambda = 1.0 - lambda;
	}

	double delta1 = DistanceChange(Previous(edge));
	double delta2 = -DistanceChange(edge);

	if( delta1 <= 0.0 )
	{
		if( delta2 <= 0.0 )
		{
			edge = Reverse(Adjacent(Previous(edge)));
			lambda = 0.0;
		} else
		{
			edge = Adjacent(Next(edge));
			lambda = 0.0;
		}
	} else
	{
		if( delta1 * lambda + delta2 <= 0.0 )
		{
			edge = Reverse(Adjacent(Previous(edge)));
			lambda = - delta1 / delta2 * lambda;
		}else
		{
			double lambda1 = (lambda * delta1 + delta2)/(delta1 + delta2);
			double delta3 = - DistanceChange(Next(edge));
			double delta4 = DistanceChange(Next(Next(edge)));
			if( (1.0-lambda1)*delta3 - lambda1*delta4 < 0 )
			{
				edge = Adjacent(Next(edge));
				lambda = (1.0 + delta3/delta4)*(1.0-lambda1);
			} else
			{
				edge = Reverse(Adjacent(Next(Next(edge))));
				lambda = (1.0 + delta4/delta3)*lambda1;
			}
		}
	}
}


Edge EikonalGeodesicRetriever::Adjacent(const Edge & edge) const
{
	Edge adj;
	adj.x = edge.x;
	adj.belongsToLeftFace = !edge.belongsToLeftFace;
	return adj;
}

Edge EikonalGeodesicRetriever::Reverse(const Edge & edge) const
{
	Edge rev;
	rev.x[0] = edge.x[1];
	rev.x[1] = edge.x[0];
	rev.belongsToLeftFace = !edge.belongsToLeftFace;
	return rev;
}

Edge EikonalGeodesicRetriever::Next(const Edge & edge) const
{
	Edge next;
	next.belongsToLeftFace = edge.belongsToLeftFace;
	next.x[0] = edge.x[1];
	next.x[1][0] = ProperMod(edge.x[1][0] + (edge.belongsToLeftFace ? 1 : -1 ) * (edge.x[0][1] - edge.x[1][1]),width_);
	next.x[1][1] = ProperMod(edge.x[1][1] + (edge.belongsToLeftFace ? 1 : -1 ) * (edge.x[1][0] - edge.x[0][0]),width_);
	return next;
}

Edge EikonalGeodesicRetriever::Previous(const Edge & edge) const
{
	Edge prev;
	prev.belongsToLeftFace = edge.belongsToLeftFace;
	prev.x[1] = edge.x[0];
	prev.x[0][0] = ProperMod(edge.x[0][0] + (edge.belongsToLeftFace ? 1 : -1 ) * ( edge.x[0][1] - edge.x[1][1] ),width_);
	prev.x[0][1] = ProperMod(edge.x[0][1] + (edge.belongsToLeftFace ? 1 : -1 ) * ( edge.x[1][0] - edge.x[0][0] ),width_);
	return prev;
}

double EikonalGeodesicRetriever::DistanceChange(const Edge & edge) const
{
	return distance_[edge.x[1][0]][edge.x[1][1]] - distance_[edge.x[0][0]][edge.x[0][1]];
}

TwoPointFunction::TwoPointFunction(VolumeMeasure * const measure, DISTANCE_METHOD method, double maxdist, int bins) :
	measure_(measure), 
	width_(measure->getWidth()),
	distance_(measure->getWidth(),std::vector<double>(measure->getWidth())),
	max_dist_(maxdist),
	bins_(bins),
	measurements_(0),
	method_(method)
{
	dist_histogram_.resize(bins_,0.0);
}

double TwoPointFunction::getDistance(int x, int y) const
{
	return distance_[x][y];
}

void TwoPointFunction::Measure(RNGenerator & rng, int n)
{
	for(int i=0;i<n;i++)
	{
		std::pair<int,int> x = measure_->RandomInMeasure(rng);
		if( method_ == EIKONAL )
		{
			EikonalDistance(x.first,x.second,true);
		} else
		{
			DijkstraDistance(x.first,x.second,true);
		}
		measurements_++;
	}
}

void TwoPointFunction::DijkstraDistance(int startX, int startY,bool measurement)
{
	for(int i=0;i<width_;i++)
	{
		for(int j=0;j<width_;j++)
		{
			distance_[i][j] = 1000.0;
		}
	}
	std::priority_queue<Node> pqueue;
	Node startNode;
	startNode.distance = 0.0;
	startNode.x[0] = startX;
	startNode.x[1] = startY;
	pqueue.push( startNode );
	distance_[startX][startY] = 0.0;

	while( !pqueue.empty() )
	{
		Node node = pqueue.top();
		pqueue.pop();

		if( node.distance > distance_[node.x[0]][node.x[1]] )
			continue;

		if( measurement && (node.x[0] != startX || node.x[1] != startY) )
		{
			if( node.distance < max_dist_ )
			{
				int bin = static_cast<int>((node.distance/max_dist_)*bins_);
				dist_histogram_[bin] += measure_->getMeasure(node.x[0],node.x[1]);
			} else
			{
				break;
			}
		}

		for(int i=0;i<4;i++)
		{
			Node nextNode(node);
			if( i==0 )
				nextNode.x[0] = (nextNode.x[0]+1)%width_;
			else if( i==1 )
				nextNode.x[0] = (nextNode.x[0]+width_-1)%width_;
			else if( i==2 )
				nextNode.x[1] = (nextNode.x[1]+1)%width_;
			else
				nextNode.x[1] = (nextNode.x[1]+width_-1)%width_;

			nextNode.distance += 0.5 * ( measure_->getSqrtMeasure(node.x[0],node.x[1]) + measure_->getSqrtMeasure(nextNode.x[0],nextNode.x[1]) );

			if( nextNode.distance < distance_[nextNode.x[0]][nextNode.x[1]] )
			{
				distance_[nextNode.x[0]][nextNode.x[1]] = nextNode.distance;
				pqueue.push( nextNode );
			}
		}
	}
}

void TwoPointFunction::EikonalDistance(int startX, int startY,bool measurement)
{
	for(int i=0;i<width_;i++)
	{
		for(int j=0;j<width_;j++)
		{
			distance_[i][j] = 100000.0;
		}
	}

	std::priority_queue<Node> pqueue;
	Node startNode;
	startNode.distance = 0.0;
	startNode.x[0] = startX;
	startNode.x[1] = startY;
	pqueue.push( startNode );
	distance_[startX][startY] = 0.0;
	double current_distance_ = 0.0;

	while( !pqueue.empty() )
	{
		Node node = pqueue.top();
		pqueue.pop();

		if( node.distance > distance_[node.x[0]][node.x[1]] )
			continue;

		current_distance_ = node.distance;

		if( measurement && (node.x[0] != startX || node.x[1] != startY) )
		{
			if( node.distance < max_dist_ )
			{
				int bin = static_cast<int>((node.distance/max_dist_)*bins_);
				dist_histogram_[bin] += measure_->getMeasure(node.x[0],node.x[1]);
			} else
			{
				break;
			}
		}

		for(int i=0;i<4;i++)
		{
			Node nextNode(node);
			std::pair<int,int> nextx = Neighbour(std::pair<int,int>(node.x[0],node.x[1]),i);
			nextNode.x[0] = nextx.first;
			nextNode.x[1] = nextx.second;

			if( distance_[nextNode.x[0]][nextNode.x[1]] > current_distance_ )
			{
				// its distance is not fixed yet
				std::pair<int,int> nextx2 = Neighbour(nextx,i);
				double a = std::min( node.distance, distance_[nextx2.first][nextx2.second] );
				std::pair<int,int> nextx3 = Neighbour(nextx,2*(1-i/2));
				std::pair<int,int> nextx4 = Neighbour(nextx,2*(1-i/2)+1);
				double b = std::min( distance_[nextx3.first][nextx3.second], distance_[nextx4.first][nextx4.second] );
				
				double localDistance = measure_->getSqrtMeasure(nextNode.x[0],nextNode.x[1]);
				if( a > 10000.0 || b > 10000.0 )
				{
					nextNode.distance = std::min(a,b) + localDistance;
				} else
				{
					double discriminant = 2 * localDistance * localDistance - (a-b)*(a-b);
					if( discriminant >= 0.0 )
					{
						nextNode.distance = 0.5*(a+b+std::sqrt(discriminant));
					} else
					{
						nextNode.distance = std::min(a,b) + localDistance;
					}
				}


				if( nextNode.distance < distance_[nextNode.x[0]][nextNode.x[1]] )
				{
					distance_[nextNode.x[0]][nextNode.x[1]] = nextNode.distance;
					pqueue.push( nextNode );
				}
			}
		}
	}
}

std::pair<int,int> TwoPointFunction::Neighbour(const std::pair<int,int> & x, int dir)
{
	return std::pair<int,int>( (dir > 1 ? x.first : (dir == 0 ? (x.first+1)%width_ : (x.first+width_-1)%width_ ) ),
		(dir <= 1 ? x.second : (dir == 2 ? (x.second+1)%width_ : (x.second+width_-1)%width_ ) ) );
}


std::string TwoPointFunction::Output()
{
	std::ostringstream os;
	os << std::fixed << "{width -> " << width_;
	os << ", method -> \"" << (method_==EIKONAL ? "Eikonal" : "Dijkstra" ) << "\"";
	os << ", maxdist -> " << max_dist_ << ", measurements -> " << measurements_;
	os << ", disthistogram -> ";
	PrintToStream(os,dist_histogram_.begin(),dist_histogram_.end());
	os << "}";
	return os.str();
}


void TwoPointFunction::SaveDistance(std::string filename) const
{
	std::ofstream file(filename.c_str());
	file << width_ << "\n";
	for(int i=0;i<width_;i++)
	{
		for(int j=0;j<width_;j++)
		{
			file << measure_->getMeasure(i,j) << "\n";
		}
	}
	for(int i=0;i<width_;i++)
	{
		for(int j=0;j<width_;j++)
		{
			file << distance_[i][j] << "\n";
		}
	}
}

CorrelationFunction::CorrelationFunction(VolumeMeasure * const measure) :
	measure_(measure),
	width_(measure->getWidth()),
	correlation_measurements_(0)
{

}

void CorrelationFunction::MeasureCorrelation(RNGenerator & rng,int n)
{
	if( static_cast<int>(correlation_.size()) < width_/2 )
	{
		correlation_.resize(width_/2,0.0);
	}
	for(int i=0;i<n;i++)
	{
		std::pair<int,int> point = measure_->RandomInMeasure(rng);
		for(int x=-(width_/2);x<width_-(width_/2);x++)
		{
			for(int y=-(width_/2);y<width_-(width_/2);y++)
			{
				int intdistance = static_cast<int>(0.5+std::sqrt(static_cast<double>(x*x + y*y)));
				if( intdistance < width_/2 )
				{
					correlation_[intdistance] += measure_->getMeasure( ProperMod(point.first + x,width_), ProperMod(point.second + y,width_) );
				}
			}
		}
		correlation_measurements_++;
	}
}

std::string CorrelationFunction::Output()
{
	std::vector<int> perbin(width_/2,0);
	for(int x=-(width_/2);x<width_-(width_/2);x++)
	{
		for(int y=-(width_/2);y<width_-(width_/2);y++)
		{
			int intdistance = static_cast<int>(0.5+std::sqrt(static_cast<double>(x*x + y*y)));
			if( intdistance < width_/2 )
			{
				perbin[intdistance]++;
			}
		}
	}

	std::ostringstream os;
	os << std::fixed << "{width -> " << width_ << ", correlationmeasurements -> " << correlation_measurements_ << ", correlation -> ";
	PrintToStream(os,correlation_.begin(),correlation_.end());
	os << ", pointsperbin -> ";
	PrintToStream(os,perbin.begin(),perbin.end());
	os << "}";
	return os.str();
}