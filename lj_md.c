#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>
#include <math.h>

// Useful references for the LJ model we're using:
// Classical Mechanics: https://www.ks.uiuc.edu/Training/Workshop/SanFrancisco/lectures/Wednesday-ForceFields.pdf (page 4)
// LJ Potential: https://en.wikipedia.org/wiki/Lennard-Jones_potential

double BOX_LOWER = 0.0;
double BOX_UPPER = 20.0;
double TIMESTEPS = 100000;
double DT = 0.001;
double EPSILON = .1;
double RADIUS = 10.0;

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
  double vx, vy, vz;
  double fx, fy, fz;
  double mass;
};

typedef struct atomStruct Atom;

struct moleculeStruct {
  int count;
  Atom** atoms;
};

typedef struct moleculeStruct Molecule;

void bound(double* value, double upper, double lower) {
  if (upper <= lower) {
    printf("Cannot have upper <= lower on a bound");
    exit(1);
  }
  double size = upper - lower;
  double res = fmod(*value, size);
  if (res >= 0)
    *value = lower + res;
  else
    *value = upper + res;
}

void update_atom_positions(Atom* a) {
  a->x += a->vx * DT;
  a->y += a->vy * DT;
  a->z += a->vz * DT;
  bound(&a->x, BOX_UPPER, BOX_LOWER);
  bound(&a->y, BOX_UPPER, BOX_LOWER);
  bound(&a->z, BOX_UPPER, BOX_LOWER);
}

void update_atom_velocities(Atom* a) {  
  a->vx += a->fx / a->mass * DT;
  a->vy += a->fy / a->mass * DT;
  a->vz += a->fz / a->mass * DT;
}

Atom* new_Atom(char element, double x, double y, double z, double mass) {
  Atom* a = malloc(sizeof(Atom));
  a->element = element;
  a->x = x;
  a->y = y;
  a->z = z;
  a->vx = 0.0;
  a->vy = 0.0;
  a->vz = 0.0;
  a->fx = 0.0;
  a->fy = 0.0;
  a->fz = 0.0;
  a->mass = mass;
  update_atom_velocities(a);
  update_atom_positions(a);
  return a;
}

Molecule* new_Molecule(int count, Atom** atoms) {
  Molecule* m = malloc(sizeof(Molecule));
  m->count = count;
  m->atoms = atoms;
  return m;
}

void update_atom_forces(Atom* a, double* forces) {
  a->fx = forces[0];
  a->fy = forces[1];
  a->fz = forces[2];
}

double atomic_distance_squared(Atom* a1, Atom* a2) {
  double xdist = a1->x - a2->x;
  double ydist = a1->y - a2->y;
  double zdist = a1->z - a2->z;
  return xdist * xdist + ydist * ydist + zdist * zdist;
}

double sigma(Atom* a1, Atom* a2) {
  return 1.5; // Assume all Hydrogen atoms
}

// Derivative of LJ potential by R in the direction of dim
// Assumes radius has been squared
double LJ_dR(double r, double dim, double e, double s) {
  double ft = 4 * e;
  double st = 6 * dim * pow(s, 6.0) / pow(r, 4.0);
  double tt = 12 * dim * pow(s, 12.0) / pow(r, 7.0);
  return ft * (st - tt);
}

// Updates the forces in molecule mol
// Assumes radius has already been squared
void update_forces(Molecule* mol, double radius) {
  double forces[3];
  Atom* copies[7];
  for (int index = 0; index < mol->count; index++) {
    forces[0] = 0.0; forces[1] = 0.0; forces[2] = 0.0;
    Atom* a = mol->atoms[index];
    double size = BOX_UPPER - BOX_LOWER;
    for (int i = 0; i < mol->count; i++) {
      if (i == index)
        continue;
      Atom* a2 = mol->atoms[i];
      Atom* copy0 = a2;
      Atom* copy1 = new_Atom(a2->element, a2->x + size, a2->y, a2->z, a2->mass);
      Atom* copy2 = new_Atom(a2->element, a2->x - size, a2->y, a2->z, a2->mass);
      Atom* copy3 = new_Atom(a2->element, a2->x, a2->y + size, a2->z, a2->mass);
      Atom* copy4 = new_Atom(a2->element, a2->x, a2->y - size, a2->z, a2->mass);
      Atom* copy5 = new_Atom(a2->element, a2->x, a2->y, a2->z + size, a2->mass);
      Atom* copy6 = new_Atom(a2->element, a2->x, a2->y, a2->z - size, a2->mass);
     
      for (int i = 0; i < 7; i++)
        copies[i] = NULL;
      if (atomic_distance_squared(a, copy0) <= radius)
        copies[0] = copy0;
      if (atomic_distance_squared(a, copy1) <= radius)
        copies[1] = copy1;
      if (atomic_distance_squared(a, copy2) <= radius)
        copies[2] = copy2;
      if (atomic_distance_squared(a, copy3) <= radius)
        copies[3] = copy3;
      if (atomic_distance_squared(a, copy4) <= radius)
        copies[4] = copy4;
      if (atomic_distance_squared(a, copy5) <= radius)
        copies[5] = copy5;
      if (atomic_distance_squared(a, copy6) <= radius)
        copies[6] = copy6;
      double s = sigma(a, mol->atoms[i]);
      for (int j = 0; j < 7; j++) {
        // != NULL
        if (copies[j]) {
          double r = atomic_distance_squared(a, copies[j]);
          forces[0] -= LJ_dR(r, a->x - copies[j]->x, EPSILON, s);
          forces[1] -= LJ_dR(r, a->y - copies[j]->y, EPSILON, s);
          forces[2] -= LJ_dR(r, a->z - copies[j]->z, EPSILON, s);
        }
      }
    }
    update_atom_forces(a, forces);
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
    atoms[i] = new_Atom(element, x, y, z, 1.0); // Assume all Hydrogen atoms
  }
  fclose(fp);
  return new_Molecule(count, atoms);
}

void write_xyz(const char* filename, Molecule* mol) {
  FILE *fp;
  fp = fopen(filename, "w");
  fprintf(fp, "%d\n\n", mol->count);
  Atom* atom;
  for (int i = 0; i < mol->count; i++) {
    atom = mol->atoms[i];
    fprintf(fp, "%c %f %f %f\n", atom->element, atom->x, atom->y, atom->z);
  }
  fclose(fp);
}

void write_xyz_format(const char* filename, Molecule* molecules, int count) {
  char f[500];
  char istr[50];
  strcpy(f, filename);
  sprintf(istr, "%06d", count);
  strcat(f, istr);
  strcat(f, ".xyz");
  write_xyz(f, molecules);
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
  srand ( time ( NULL));
  double forces[3];
  int write_frequency = (int) (.1 / DT);
  for (int i = 0; i < TIMESTEPS; i++) {
    if (i % (write_frequency * 100) == 0)
      printf("Timestep %d\n", i);
    if (i % write_frequency == 0)
      write_xyz_format(outstr, mol, i / write_frequency);
    update_forces(mol, RADIUS * RADIUS);
    for (int j = 0; j < mol->count; j++) {
      update_atom_velocities(mol->atoms[j]);
      update_atom_positions(mol->atoms[j]);
    }
  }
  
  return 0;
}