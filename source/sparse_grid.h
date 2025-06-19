/******************************************************************************
 * Support for a sparse grid with a random offset on each timeslice 
 * ***************************************************************************/

#ifndef SPARSE_GRID_H
#define SPARSE_GRID_H
#include "array.h"
#include "layout_info.h"
#include <random>
#include <cstdint>

namespace LaphEnv {

	class RandomSparseGrid {
		int grid_spacing, n_grid_points;
		bool is_random; 
		uint32_t seed;
		std::vector<std::vector<int>> offsets; 

		struct SparseOffset { 
			int local_offset; 
			int global_offset;

			SparseOffset(int _local_offset, int _global_offset) : 
				local_offset(_local_offset), global_offset(_global_offset) {} 
		}; 
		
		void check_grid_spacing() {
			if (grid_spacing<0) {
				throw std::logic_error("invalid grid spacing");
			}
			if (((grid_spacing%LayoutInfo::getLattExtents()[0])!=0)||
					((grid_spacing%LayoutInfo::getLattExtents()[1])!=0)||
					((grid_spacing%LayoutInfo::getLattExtents()[2])!=0)) {
				throw std::logic_error("grid spacing must divide local spatial extents");
			}
		}	

		void set_number_of_grid_points() {
			int Lx = LayoutInfo::getLattExtents()[0];
			int Ly = LayoutInfo::getLattExtents()[1];
			int Lz = LayoutInfo::getLattExtents()[2];

			n_grid_points = Lx*Ly*Lz/(grid_spacing*grid_spacing*grid_spacing); 
		}

		void generate_offsets() {
  		offsets.clear();     
			int Lt = LayoutInfo::getLattExtents()[3];
      if (is_random) {
        std::mt19937 gen(seed);
        std::uniform_int_distribution<int> dist(0,grid_spacing-1);
        for (int t=0;t<Lt;t++) {
          std::vector<int> coords;
          coords.push_back(dist(gen));
          coords.push_back(dist(gen));
          coords.push_back(dist(gen));
          offsets.push_back(coords);
        }
      } else {
        std::vector<int> coords = {0, 0, 0};
        for (int t=0;t<Lt;t++) {
          offsets.push_back(coords);
        }
      }
    }


		public:

		RandomSparseGrid(int _grid_spacing, const uint32_t& _seed=0) : 
			grid_spacing(_grid_spacing), seed(_seed), is_random(true) {	
        if (seed==0) 
					is_random=false;
				check_grid_spacing(); 
				generate_offsets();
				set_number_of_grid_points();
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
			check_grid_spacing();
			generate_offsets();
			set_number_of_grid_points();
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

		//Return the offsets corresponding to the sparse grid points on the local lattice 
		std::vector<SparseOffset> getLocalGridPoints(int t) const {
			if ((t<0)||(t>=LayoutInfo::getLattExtents()[3]))
				throw std::logic_error("invalid time Extent in getLocalGridPoints");

			std::vector<SparseOffset> ret;
			std::vector<int> coords = {0,0,0,t};
		  int Lx = LayoutInfo::getLattExtents()[0];	
		  int Ly = LayoutInfo::getLattExtents()[1];	
		  int Lz = LayoutInfo::getLattExtents()[2];

			int ctr; 
			for (int iZ=0; iZ<Lz/grid_spacing; iZ++) {
				coords[2]=(iZ*grid_spacing+offsets[t][2])%Lz;
				for (int iY=0; iY<Ly/grid_spacing; iY++) {
					coords[1]=(iY*grid_spacing+offsets[t][1])%Ly;
					for (int iX=0; iX<Lx/grid_spacing; iX++) {
						coords[0]=(iX*grid_spacing+offsets[t][0])%Lx;
						int rank,rank_site_linear_index; 
						LayoutInfo::getCommInfoFromLatticeCoords(
								coords, rank, rank_site_linear_index);
						if (rank==LayoutInfo::getMyRank())
							ret.push_back(SparseOffset(rank_site_linear_index,ctr));
						ctr++;
					}
				}
			}
			return ret; 
		}

		int getGridSpacing() const { return grid_spacing; }
		int getNGridPoints() const { return n_grid_points; }
    bool isRandom() const { return is_random; }

		private:
		// prevent copying
		RandomSparseGrid(const RandomSparseGrid&);
		RandomSparseGrid& operator=(const RandomSparseGrid&);
	};

	// ******************************************************************
}
#endif
