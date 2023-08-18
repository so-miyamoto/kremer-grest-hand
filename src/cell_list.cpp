#include "cell_list.h"


CellList::CellList()
{
  _initialized_boxsize = false;
  _cutted_box = false;
};

void 
CellList::set_boxsize(double Lx,double Ly, double Lz)
{
  _initialized_boxsize = true;
  _Lx = Lx;
  _Ly = Ly;
  _Lz = Lz;
};

void 
CellList::cut_box(double minimum_cell_length){
  _cutted_box = true;
  if( !_initialized_boxsize ){
    throw std::logic_error("CellList::cannot determine cell sizes before initialize box\n");
  }
  if( _Lx < 3*minimum_cell_length ){
    throw std::runtime_error("CellList::too small to cut Lx\n");
  }
  if( _Ly < 3*minimum_cell_length ){
    throw std::runtime_error("CellList::too small to cut Ly\n");
  }
  if( _Lz < 3*minimum_cell_length ){
    throw std::runtime_error("CellList::too small to cut Lz\n");
  }

  _num_cells_x = int(_Lx/minimum_cell_length) - 1;
  _num_cells_y = int(_Ly/minimum_cell_length) - 1;
  _num_cells_z = int(_Lz/minimum_cell_length) - 1;

  _cellLx = _Lx/_num_cells_x;
  _cellLy = _Ly/_num_cells_y;
  _cellLz = _Lz/_num_cells_z;
  
  _num_cells = _num_cells_x*_num_cells_y*_num_cells_z;
  _cells_list.resize(_num_cells);

}

void
CellList::register_addresses(std::vector<Vecd>& pos){
  if( !_cutted_box ){
    throw std::logic_error("CellList::cannot distribute indexes to cells before init\n");
  }

  for (int i = 0; i < _num_cells; i++)
    _cells_list[i].clear();
  
  const int N = pos.size();
  for (int i = 0; i < N; i++)
  {
    auto& p = pos[i];
    int ix = p[0] / _cellLx;
    int iy = p[1] / _cellLy;
    int iz = p[2] / _cellLz;
    for(int nx = -1; nx < 2; nx++)
      for(int ny = -1; ny < 2; ny++)
        for(int nz = -1; nz < 2; nz++)
          _cells_list[_idx(ix+nx,iy+ny,iz+nz)].emplace_back(i);
  }
  return;
}

std::vector<int>*
CellList::get_neighbors(const Vecd p){
  int ix = p[0] / _cellLx;
  int iy = p[1] / _cellLy;
  int iz = p[2] / _cellLz;
  return &(_cells_list[_idx(ix,iy,iz)]);
}



int 
CellList::_idx(int ix, int iy, int iz) const {
  ix = (ix+_num_cells_x)%_num_cells_x;
  iy = (iy+_num_cells_y)%_num_cells_y;
  iz = (iz+_num_cells_z)%_num_cells_z;    
  return ix*_num_cells_y*_num_cells_z + iy*_num_cells_z + iz;
}



/** test

Vecd boundary_condition(Vecd p1,Vecd p2,double Lx,double Ly,double Lz){
  constexpr int x=0,y=1,z=2;
  Vecd q = p1-p2;
  if (q[x] < -0.5*Lx) q[x] += Lx; // left x
  if (q[y] < -0.5*Ly) q[y] += Ly; // left y
  if (q[z] < -0.5*Lz) q[z] += Lz; // left z

  if (q[x] >= 0.5*Lx) q[x] -= Lx; // right x
  if (q[y] >= 0.5*Ly) q[y] -= Ly; // right y
  if (q[z] >= 0.5*Lz) q[z] -= Lz; // right z  
  return q;
}

int main(){

  try{
  CellList cell_list;

  cell_list.set_boxsize(10.0,20.0,30.0);
  cell_list.cut_box(2.0);

  std::vector<Vecd> pos(10);
  for (int i = 0; i < 10; i++)
  {
    double ix = i;
    pos[i] = Vecd{ix,ix,ix};
    std::cout<<pos[i]<<std::endl;
  }

  cell_list.register_addresses(pos);
  for (size_t i = 0; i < pos.size(); i++)
  {
    std::vector<int>* lists = cell_list.get_neighbors(pos[i]);
    std::cout<<"i="<<i<<"\n";
    for (size_t j = 0; j < lists->size(); j++)
    {
      std::cout<<(*lists)[j]<<", dist="<<norm(boundary_condition(pos[i],pos[(*lists)[j]],10.0,20.0,30.0))<<"\n";
    }
    std::cout<<std::endl;
    
  }
  }catch(std::runtime_error& e){
    std::cout<<e.what()<<std::endl;
  }catch(std::logic_error& e){
    std::cout<<e.what()<<std::endl;
  }
  return 0;  

}


*/