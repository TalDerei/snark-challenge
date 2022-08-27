// This directory contains a reference CPU implementation of Field Arithmetic using libff (C++ library for Finite Fields and Elliptic Curves)

#include <cstdio>
#include <vector>

#include <libff/algebra/curves/mnt753/mnt4753/mnt4753_pp.hpp>
#include <libff/algebra/curves/mnt753/mnt6753/mnt6753_pp.hpp>
#include <libff/algebra/curves/mnt753/mnt4753/mnt4753_init.hpp>

using namespace libff;
using namespace std;

Fq<mnt4753_pp> read_mnt4_fq(FILE* input) {
  // Define a finite field with q elements : mnt4753_pp
  Fq<mnt4753_pp> x;
  Fq<mnt4753_pp> y;
  // Defining pointer to array of 12 64-bit limbs, mnt4753_q_limbs is 12L representing the size 
  // of each element to be read * long int, and input is a file pointer to the data stream.
  // This is represented as a little-endian length 12 array of 64-bit integers
  fread((void *) x.mont_repr.data, libff::mnt4753_q_limbs * sizeof(mp_size_t), 1, input);
  return x;
}

void write_mnt4_fq(FILE* output, Fq<mnt4753_pp> x) {
  fwrite((void *) x.mont_repr.data, libff::mnt4753_q_limbs * sizeof(mp_size_t), 1, output);
}

Fq<mnt6753_pp> read_mnt6_fq(FILE* input) {
  // bigint<mnt4753_q_limbs> n;
  Fq<mnt6753_pp> x;
  fread((void *) x.mont_repr.data, libff::mnt6753_q_limbs * sizeof(mp_size_t), 1, input);
  return x;
}

void write_mnt6_fq(FILE* output, Fq<mnt6753_pp> x) {
  fwrite((void *) x.mont_repr.data, libff::mnt6753_q_limbs * sizeof(mp_size_t), 1, output);
}

Fq<mnt4753_pp> read_mnt4_fq_numeral(FILE* input) {
  // bigint<mnt4753_q_limbs> n;
  Fq<mnt4753_pp> x;
  fread((void *) x.mont_repr.data, libff::mnt4753_q_limbs * sizeof(mp_size_t), 1, input);
  auto b = Fq<mnt4753_pp>(x.mont_repr);
  return b;
}

void write_mnt4_fq_numeral(FILE* output, Fq<mnt4753_pp> x) {
  auto out_numeral = x.as_bigint();
  fwrite((void *) out_numeral.data, libff::mnt4753_q_limbs * sizeof(mp_size_t), 1, output);
}

Fq<mnt6753_pp> read_mnt6_fq_numeral(FILE* input) {
  // bigint<mnt4753_q_limbs> n;
  Fq<mnt6753_pp> x;
  fread((void *) x.mont_repr.data, libff::mnt6753_q_limbs * sizeof(mp_size_t), 1, input);
  auto b = Fq<mnt6753_pp>(x.mont_repr);
  return b;
}

void write_mnt6_fq_numeral(FILE* output, Fq<mnt6753_pp> x) {
  auto out_numeral = x.as_bigint();
  fwrite((void *) out_numeral.data, libff::mnt6753_q_limbs * sizeof(mp_size_t), 1, output);
}

void print_array(uint8_t* a) {
  for (int j = 0; j < 96; j++) {
    printf("%x ", ((uint8_t*)(a))[j]);
  }
  printf("\n");
}

// The actual code for doing Fq multiplication lives in libff/algebra/fields/fp.tcc
int main(int argc, char *argv[])
{
    // argv should be
    // { "main", "compute" or "compute-numeral", inputs, outputs }

    mnt4753_pp::init_public_params();
    mnt6753_pp::init_public_params();

    size_t n;

    auto is_numeral = strcmp(argv[1], "compute-numeral") == 0;
    auto inputs = fopen(argv[2], "r");
    auto outputs = fopen(argv[3], "w");

    auto read_mnt4 = read_mnt4_fq;
    auto read_mnt6 = read_mnt6_fq;
    auto write_mnt4 = write_mnt4_fq;
    auto write_mnt6 = write_mnt6_fq;
    if (is_numeral) {
      read_mnt4 = read_mnt4_fq_numeral;
      read_mnt6 = read_mnt6_fq_numeral;
      write_mnt4 = write_mnt4_fq_numeral;
      write_mnt6 = write_mnt6_fq_numeral;

    }

    // The while loop increments n from 1024,2048,4096,8192 creating arrays for both
    // the MNT4-753 curve and MNT6-753 curve, doing a multiplication for each respective curve 
    // and adding the result to an output file. 
    while (true) {
      size_t elts_read = fread((void *) &n, sizeof(size_t), 1, inputs);
      // cout << "size of elts_read is: " << elts_read << endl;
      if (elts_read == 0) { break; }

      // x0, x1, y0, y1 : Array(𝔽MNT4753.q, n)
      std::vector<Fq<mnt4753_pp>> x0;
      cout << "size of n is: " << n << endl;

      for (size_t i = 0; i < n; ++i) {
        x0.emplace_back(read_mnt4(inputs));
      }

      std::vector<Fq<mnt4753_pp>> x1;
      for (size_t i = 0; i < n; ++i) {
        x1.emplace_back(read_mnt4(inputs));
      }

      std::vector<Fq<mnt6753_pp>> y0;
      for (size_t i = 0; i < n; ++i) {
        y0.emplace_back(read_mnt6(inputs));
      }
      std::vector<Fq<mnt6753_pp>> y1;
      for (size_t i = 0; i < n; ++i) {
        y1.emplace_back(read_mnt6(inputs));
      }

      // out_x, out_y : Array(𝔽MNT4753.q, n)
      // output out_x should have out_x[i] = x0[i] * x1[i]. where * is multiplication in the field 𝔽MNT4753.q.
      // output out_y should have out_y[i] = y0[i] * y1[i]. where * is multiplication in the field 𝔽MNT6753.q.
      for (size_t i = 0; i < n; ++i) {
        write_mnt4(outputs, x0[i] * x1[i]);
        // cout << "x0 is " << endl;
        // x0[i].print();
        // cout << "x1 is " << endl;
        // x1[i].print();
        // cout << "multiplication is " << endl;
        // auto sum  = x0[i] * x1[i];
        // sum.print();
        // exit(0);
      }

      for (size_t i = 0; i < n; ++i) {
        write_mnt6(outputs, y0[i] * y1[i]);
      }
    }
    fclose(outputs);

    return 0;
}
