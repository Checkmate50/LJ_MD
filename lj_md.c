#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

struct atomStruct {
  char element;
  double x, y, z;
};

struct moleculeStruct {
  int size;
  struct Atom* atoms;
};

typedef struct atomStruct Atom;
typedef struct moleculeStruct Molecule;

Atom* new_Atom(char element, double x, double y, double z) {
  Atom* a = malloc(sizeof(Atom));
  a->element = element;
  a->x = x;
  a->y = y;
  a->z = z;
  return a;
}

Molecule* new_Molecule(int size, struct Atom* atoms) {
  Molecule* m = malloc(sizeof(Molecule));
  m->size = size;
  m->atoms = atoms;
  return m;
}

void slice_str(const char * str, char * buffer, size_t start, size_t end)
{
  size_t j = 0;
  for ( size_t i = start; i <= end; ++i ) {
    buffer[j++] = str[i];
  }
  buffer[j] = 0;
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
  for (int i = 0; i < count; i++) {
    fgets(buff, 255, fp);
    char element = strtok(buff, deliminator)[0];
    double x = strtod(strtok(NULL, deliminator), &END);
    double y = strtod(strtok(NULL, deliminator), &END);
    double z = strtod(strtok(NULL, deliminator), &END);
    printf("%c %f %f %f\n", element, x, y, z);
  }
  fclose(fp);
}

void write_xyz(const char* filename, Molecule* mol) {
  FILE *fp;
  fp = fopen(filename, "w");
  fclose(fp);
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
  write_xyz(outstr, mol);
  
  return 0;
}