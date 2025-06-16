/******************************************************************************
 * Support for a sparse grid with a random offset on each timeslice 
 * ***************************************************************************/

#ifndef SPARSE_GRID_H
#define SPARSE_GRID_H
#include "array.h"
#include <random>
#include <cstdint>

namespace LaphEnv {

class RandomSparseGrid {
	int grid_spacing;
	bool is_random; 
  uint32_t seed;
	public:

	RandomSparseGrid(int _grid_spacing, const uint32_t& _seed) : grid_spacing(grid_spacing), seed(_seed), is_random(true) {
		if (grid_spacing<0) { 
			throw std::logic_error("invalid grid spacing"); 
		}
	}

		RandomSparseGrid(int _grid_spacing) : grid_spacing(grid_spacing), is_random(false), seed(0) {
		if (grid_spacing<0) { 
			throw std::logic_error("invalid grid spacing"); 
		}
	}


	RandomSparseGrid(const XMLHandler& xmlin) : is_random(false), seed(0) { 
		XMLHandler xml_in(xmlin);
		xml_tag_assert(xml_in,"RandomSparseGridInfo","RandomSparseGrid");
		XMLHandler xmlr(xml_in, "RandomSparseGridInfo");
		if (xml_tag_count(xml_in,"RandomSeed")>0) {
			is_random=true;
			xmlread(xml_in,"RandomSeed", seed, "RandomSparseGrid");
		}
		xmlread(xml_in,"GridSpacing", grid_spacing, "RandomSparseGrid");
		if (grid_spacing<0){ 
			xmlreadfail(xml_in,"RandomSparseGrid",
					"invalid grid_spacing in RandomSparseGrid");}
	}

	std::string output(int indent = 0) const { 
		XMLHandler xmlout;
		output(xmlout);
		return xmlout.output(indent);
	}

	void output(XMLHandler& xmlout) const {
		xmlout.set_root("RandomSparseGridInfo");
		xmlout.put_child("RandomSeed", make_string(seed));
		xmlout.put_child("GridSpacing", make_string(grid_spacing));
	}

	void checkEqual(const RandomSparseGrid& in) const {
		if  ((seed!=in.seed)||(grid_spacing!=in.grid_spacing)){
    std::cerr << "RandomSparseGrid checkEqual failed"<<std::endl;
    std::cerr << "LHS:"<<std::endl<<output()<<std::endl<<"RHS:"<<std::endl<<in.output()<<std::endl;
    throw(std::invalid_argument("RandomSparseGrid checkEqual failed..."));}
	}

	std::vector<std::vector<int>> generateOffsets(int Lx, int Ly, int Lz, 
			int Lt) const {
		if (Lt<0) 
			throw std::logic_error("invalid time extent"); 

		if (((Lx%grid_spacing)!=0)||((Ly%grid_spacing)!=0)||
				((Lz%grid_spacing)!=0)) 
			throw std::logic_error("incompatible spatial lattice dimensions"); 

		std::vector<std::vector<int>> ret; 

		if (is_random) {
			std::mt19937 gen(seed);
			std::uniform_int_distribution<int> distX(0,Lx-1);
			std::uniform_int_distribution<int> distY(0,Ly-1);
			std::uniform_int_distribution<int> distZ(0,Lz-1);
			for (int t=0;t<Lt;t++) {
			        std::vector<int> coords;	
				coords.push_back(distX(gen));
				coords.push_back(distY(gen));
				coords.push_back(distZ(gen));
				ret.push_back(coords);
			} 
		}	
		else {
			std::vector<int> coords; 
			coords.push_back(0);
			coords.push_back(0);
			coords.push_back(0);
			for (int t=0;t<Lt;t++) {
				ret.push_back(coords);	
			}
		}
		return ret;
	}

	int getGridSpacing() const { return grid_spacing; }

	private:
	// prevent copying
	RandomSparseGrid(const RandomSparseGrid&);
	RandomSparseGrid& operator=(const RandomSparseGrid&);
};

// ******************************************************************
}
#endif
