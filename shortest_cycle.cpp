#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

#include "liouville.h"
#include "utilities.h"

int main(int argc, char* argv[])
{
	ParameterStream param(argc,argv);

	int width =	param.Read<int>("width"); 
	double gamma = param.Read<double>("gamma");
	int method = param.Read<int>("averaging method (0=BOX_DIVISION,1=BOX_AVERAGING,2=DISK_AVERAGING,3=DISK_DIFFUSION)");
	double deltamin = param.Read<double>("delta min");
	double deltamax = param.Read<double>("delta max");
	double deltafactor = param.Read<double>("delta factor");
	bool normalize = (param.Read<int>("normalize volume (1=yes,0=no)")==1); 
	int sweepsperoutput = param.Read<int>("measurements per delta per output");
	std::string path = param.Read<std::string>("output path");
	bool output = param.UserInput();

	if( path.length() == 0 )
	{
		path = ".";
	}
	if( path.back() != '/' && path.back() != '\\' )
	{
		path += '/';
	}

	RNGenerator rng;
	GaussianField gaussian(width);
	VolumeMeasure measure(width);
	VolumeMeasure averagedmeasure(width);
	std::vector<double> deltas;
	std::vector<ShortestCycle *> cycle;
	double maxdist = 1.5 - 0.3*gamma;
	for(double delta = deltamin;delta <= 1.00001 * deltamax; delta *= deltafactor )
	{
		deltas.push_back(delta);
		cycle.push_back(new ShortestCycle(&averagedmeasure,ShortestCycle::EIKONAL,maxdist));
	}

	VolumeMeasure::AveragingMethod avmethod;
	std::string avmethodstr;
	if(method==0)
	{
		avmethod = VolumeMeasure::BOX_SUBDIVISION;
		avmethodstr = "BOX_SUBDIVISION";
	} else if(method==1)
	{
		avmethod = VolumeMeasure::BOX_AVERAGING;
		avmethodstr = "BOX_AVERAGING";
	} else if(method==2)
	{
		avmethod = VolumeMeasure::DISK_AVERAGING;
		avmethodstr = "DISK_AVERAGING";
	} else if(method==3)
	{
		avmethod = VolumeMeasure::DISK_DIFFUSION;
		avmethodstr = "DISK_DIFFUSION";
	} else
	{
		std::cout << "Method unknown.\n";
		return 1;
	}

	while(true)
	{
		for(int i=0,endi=deltas.size();i<endi;i++)
		{
			double delta = deltas[i];
		
			if(output)
				std::cout << "delta = " << delta << "\n";
			for(int sweeps=0;sweeps<sweepsperoutput;sweeps++)
			{
				if(output)
					std::cout << sweeps << ": generate - ";

				gaussian.GenerateRandom(rng);
				measure.setFromGaussianField(gaussian,gamma);
				if( normalize )
				{
					measure.NormalizeToUnitVolume();
				}
				measure.performAveraging(averagedmeasure,avmethod,delta);

				if(output)
					std::cout << "cycle - ";

				double len = cycle[i]->FindLength();

				if(output)
					std::cout << len << " - done\n";
			}
		}
		std::ostringstream os;
		os << path << "cycl-" << width << "-" << gamma << "-" << method << "-" << (normalize?"1":"0") << ".txt";
		std::ofstream file(os.str().c_str());
		file << std::fixed << "{ commandline -> \"" << param.getCommandLine() << "\"";
		file << ", width -> " << width << ", gamma -> " << gamma;
		file << ", averagingmethod -> \"" << avmethodstr << "\"";
		file << ", normalized -> " << (normalize?"True":"False");
		file << ", datasets -> {";
		for(int i=0,endi=deltas.size();i<endi;i++)
		{
			file << (i>0?",":"") << "{ delta -> " << deltas[i]*1000.0 << "/1000";
			file << ", cycledata -> " << cycle[i]->Output() << "}";
	
		}
		file << "}}\n";
		if(output)
			std::cout << "Output: " << os.str() << "\n";	
	}

	return 0;
}

