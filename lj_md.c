#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

void slice_str(const char * str, char * buffer, size_t start, size_t end)
{
  size_t j = 0;
  for ( size_t i = start; i <= end; ++i ) {
    buffer[j++] = str[i];
  }
  buffer[j] = 0;
}

struct atomStruct {
  char element;
  double x, y, z;
};

typedef struct atomStruct Atom;

struct moleculeStruct {
  int size;
  Atom** atoms;
};

typedef struct moleculeStruct Molecule;

Atom* new_Atom(char element, double x, double y, double z) {
  Atom* a = malloc(sizeof(Atom));
  a->element = element;
  a->x = x;
  a->y = y;
  a->z = z;
  return a;
}

Molecule* new_Molecule(int size, Atom** atoms) {
  Molecule* m = malloc(sizeof(Molecule));
  m->size = size;
  m->atoms = atoms;
  return m;
}

void bound(double* value, double upper, double lower) {
  size = upper - lower;
  if (*value > upper)
    *value -= size;
  else if (*value < lower)
    *value += size;
}

void calculate_new_positions(Molecule* mol) {
  srand ( time ( NULL));
  
  for (int i = 0; i < mol->count; i++) {
    mol->atoms[i]->x += ((double)rand()/RAND_MAX) - .5;
    mol->atoms[i]->y += ((double)rand()/RAND_MAX) - .5;
    mol->atoms[i]->z += ((double)rand()/RAND_MAX) - .5;
    
    bound(&mol->atoms[i]->x, 5.0, 0.0);
    bound(&mol->atoms[i]->y, 5.0, 0.0);
    bound(&mol->atoms[i]->z, 5.0, 0.0);
  }
}

bool check_xyz(const char* filename) {
  int l = strlen(filename);
  if(l < 4) {
    return false;
  }
  char res[4];
  slice_str(filename, res, l - 4, l - 1);
  return !strcmp(res, ".xyz");
}

Molecule* read_xyz(const char* filename) {
  FILE *fp;
  char buff[255];
  char* END;
  fp = fopen(filename, "r");
  fgets(buff, 255, fp);
  long count = strtol(buff, &END, 10);
  
  fgets(buff, 255, fp); // Skip blank line
  char* deliminator = " \t";
  Atom** atoms = malloc(sizeof(Atom*) * count);
  
  for (int i = 0; i < count; i++) {
    fgets(buff, 255, fp);
    // Tokenize the element name and the XYZ positions
    // Expect XYZ file format, don't bother checking this stuff
    char element = strtok(buff, deliminator)[0];
    double x = strtod(strtok(NULL, deliminator), &END);
    double y = strtod(strtok(NULL, deliminator), &END);
    double z = strtod(strtok(NULL, deliminator), &END);
    atoms[i] = new_Atom(element, x, y, z);
  }
  fclose(fp);
  return new_Molecule(count, atoms);
}

void write_xyz(const char* filename, Molecule* mol) {
  FILE *fp;
  fp = fopen(filename, "w");
  fprintf(fp, "%d\n\n", mol->size);
  Atom* atom;
  for (int i = 0; i < mol->size; i++) {
    atom = mol->atoms[i];
    fprintf(fp, "%c %f %f %f\n", atom->element, atom->x, atom->y, atom->z);
  }
  fclose(fp);
}

void write_xyz_format(const char* filename, Molecule* molecules, int count) {
  char f[500];
  char istr[50];
  strcpy(f, filename);
  sprintf(istr, "%d", count);
  strcat(f, istr);
  strcat(f, ".xyz");
  write_xyz(f, molecules[i]);
}

int main(int argc, char* argv[]) {
  if (argc < 3) {
    printf("Must give 2 arguments: input.xyz and output");
    return 1;
  }
  
  char* instr = argv[1];
  char* outstr = argv[2];
  
  if (!check_xyz(instr)) {
    printf("Bad input filename %s", instr);
    return 1;
  }
  
  if (strlen(outstr) < 1) {
    printf("Bad output filename %s", outstr);
    return 1;
  }
  
  Molecule* mol = read_xyz(instr);
  for (int i = 0; i < 5; i++) {
    write_xyz_format(outstr, mol, i);
  }
  
  return 0;
}