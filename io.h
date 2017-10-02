#ifndef IO_H
#define IO_H

#include <iostream>
#include <fstream>

using namespace std;

class io{
public:
  io(string file_name, string type);
  ~io();
  void write(string s);
  void writeFloat(float d);
  void writeInt(int i);

  bool readLine(string &s);
  bool readWord(string &s);
  bool readFloat(float &n);
  bool readInt(int &n);
  bool readBool(bool &n);

  void newLine();
  void setAsReadOnly();
  void setAsWriteOnly();
  void force();
  void close();

  static constexpr const char* type_read = "read";
  static constexpr const char* type_write = "write"; //TODO: Something better? enum?



private:
  string type_;
  string file_name_;
  ofstream out;
  ifstream in;

};

void io::force(){
  out << flush;
}

io::io(string file_name, string type){
  type_ = type;
  file_name_ = file_name;
  if(type_ == type_read){
    in.open(file_name_);
    return;
  }
  if(type_ == type_write){
    out.open(file_name_);
    return;
  }
  cerr << "Error: valid io types are read and write" << endl;
}

io::~io(){
  out.close();
}

void io::write(string s){
  if(type_ == type_read){
    cerr << "Error: Trying to write to read only io" << endl;
    return;
  }
  out << s;
}

void io::writeFloat(float d){
  if(type_ == type_read){
    cerr << "Error: Trying to write to read only io" << endl;
    return;
  }
  out << d;
}

void io::writeInt(int i){
  if(type_ == type_read){
    cerr << "Error: Trying to write to read only io" << endl;
    return;
  }
  out << i;
}

void io::newLine(){
  if(type_ == type_read){
    cerr << "Error: Trying to write to read only io" << endl;
    return;
  }
  out << "\n";
}


bool io::readLine(string &s){
  if(type_ == type_write){
    cerr << "Warning: Trying to set write only file as read onlty." << endl;
    return false;
  }else{
       return (bool)getline(in,s);
  }
}

bool io::readWord(string &s){
  if(type_ == type_write){
    cerr << "Warning: Trying to set write only file as read onlty." << endl;
    return false;
  }else{
    return (bool)(in >> s);
  }
}

bool io::readFloat(float &d){
  if(type_ == type_write){
    cerr << "Warning: Trying to set write only file as read onlty." << endl;
    return false;
  }else{
    return (bool)(in >> d);
  }
}

bool io::readInt(int &i){
  if(type_ == type_write){
    cerr << "Warning: Trying to set write only file as read onlty." << endl;
    return false;
  }else{
    return (bool)(in >> i);
  }
}

bool io::readBool(bool &b){
  if(type_ == type_write){
    cerr << "Warning: Trying to set write only file as read onlty." << endl;
    return false;
  }else{
    return (bool)(in >> b);
  }
}


void io::close(){
 in.close();
 out.close();
}


#endif
