#include <iostream>
#include <fstream>

using namespace std;

class io{
public:
  io(string file_name);
  ~io();
  void write(string s);
  void newLine();
  void setAsReadOnly();

private:
  bool read_only_;
  string file_name_;
  ofstream out;
  unsigned int buff_size_;
};


io::io(string file_name){
  file_name_ = file_name;
  out.open(file_name_);
  read_only_ = false;
}

io::~io(){
  out.close();
}

void io::write(string s){
  if(read_only_){
    cerr << "Error: Trying to write to read only io" << endl;
    return;
  }
  out << s;
}

void io::newLine(){
  if(read_only_){
    cerr << "Error: Trying to write to read only io" << endl;
    return;
  }
  out << "\n";
}


void io::setAsReadOnly(){
  read_only_ = true;
}
