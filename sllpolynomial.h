/**
 * AUTOR: RAUL NAVARRO COBOS
 * FECHA: 07/04/2026
 * EMAIL: alu0101484365@ull.edu.es
 * VERSION: 1.0
 * ASIGNATURA: Algoritmo y Estructuras de Datos
 * PRACTICA: 4
 * ESTILO: Google C++ Style Guide
 * COMENTARIOS: 
 * COMPILACION: g++ -g main_sllpolynomial.cc -o main_sllpolynomial
 * EJECUCION: ./main_sllpolynomial < data_sllpolynomial.txt
 */

#ifndef SLLPOLYNOMIAL_H_
#define SLLPOLYNOMIAL_H_

#include <iostream>
#include <math.h>  // fabs, pow

#include "pair_t.h"
#include "sll_t.h"
#include "vector_t.h"

#define EPS 1.0e-6

typedef pair_t<double> pair_double_t;  // Campo data_ de SllPolynomial
typedef sll_node_t<pair_double_t> SllPolyNode;  // Nodos de SllPolynomial

// Clase para polinomios basados en listas simples de pares
class SllPolynomial : public sll_t<pair_double_t> {
 public:
  // constructores
  SllPolynomial(void) : sll_t() {};
  SllPolynomial(const vector_t<double>&, const double = EPS);

  // destructor
  ~SllPolynomial() {};

  // E/S
  void Write(std::ostream& = std::cout) const;
  
  // operaciones
  double Eval(const double) const;
  bool IsEqual(const SllPolynomial&, const double = EPS) const;
  void Sum(const SllPolynomial&, SllPolynomial&, const double = EPS);
  
  // MODIFICACION
  double suma_coef() const;
};


bool IsNotZero(const double val, const double eps = EPS) {
  return fabs(val) > eps;
}

// FASE II
// constructor
SllPolynomial::SllPolynomial(const vector_t<double>& v, const double eps) {
  for (int i = v.get_size() - 1; i >= 0; i--) {
    if (IsNotZero(v[i], eps)) {
      push_front(new SllPolyNode(pair_double_t(v[i], i)));
    }
  }
}


// E/S
void SllPolynomial::Write(std::ostream& os) const {
  os << "[ ";
  bool first{true};
  SllPolyNode* aux{get_head()};
  while (aux != NULL) {
    int inx{aux->get_data().get_inx()};
    double val{aux->get_data().get_val()};
    if (val > 0)
      os << (!first ? " + " : "") << val;
    else
      os << (!first ? " - " : "-") << fabs(val);
    os << (inx > 1 ? " x^" : (inx == 1) ? " x" : "");
    if (inx > 1)
      os << inx;
    first = false;
    aux = aux->get_next();
  }
  os << " ]" << std::endl;
}

std::ostream& operator<<(std::ostream& os, const SllPolynomial& p) {
  p.Write(os);
  return os;
}


// Operaciones con polinomios

// FASE III
// Evaluación de un polinomio representado por lista simple
double SllPolynomial::Eval(const double x) const {
  double result{0.0};
  SllPolyNode* aux = get_head();
  while (aux != NULL) {
    double coef = aux->get_data().get_val();
    int grado = aux->get_data().get_inx();
    result += coef * pow(x, grado);
    aux = aux->get_next();
  }
  
  return result;
}

// Comparación si son iguales dos polinomios representados por listas simples
bool SllPolynomial::IsEqual(const SllPolynomial& sllpol, const double eps) const {
  bool differents = false;
  SllPolyNode* a = get_head();
  SllPolyNode* b = sllpol.get_head();
  while (a != NULL && b != NULL) {
    int grado_a = a->get_data().get_inx();
    int grado_b = b->get_data().get_inx();
    double val_a = a->get_data().get_val();
    double val_b = b->get_data().get_val();
    if (grado_a != grado_b || fabs(val_a - val_b) > eps) {
      differents = true;
      break;
    }
    a = a->get_next();
    b = b->get_next();
  }

  // Si una lista termina antes que la otra → diferentes
  if (a != NULL || b != NULL) {
    differents = true;
  }
  return !differents;
}

// FASE IV
// Generar nuevo polinomio suma del polinomio invocante mas otro polinomio
void SllPolynomial::Sum(const SllPolynomial& sllpol, SllPolynomial& sllpolsum, const double eps) {
  while (!sllpolsum.empty()) {
    delete sllpolsum.pop_front();
  }
  SllPolyNode* a = get_head();
  SllPolyNode* b = sllpol.get_head();
  SllPolyNode* tail = NULL;
  while (a != NULL && b != NULL) {
    int grado_a = a->get_data().get_inx();
    int grado_b = b->get_data().get_inx();
    double val_a = a->get_data().get_val();
    double val_b = b->get_data().get_val();

    SllPolyNode* nuevo = NULL;

    if (grado_a == grado_b) {
      double suma = val_a + val_b;
      if (IsNotZero(suma, eps)) {
        nuevo = new SllPolyNode(pair_double_t(suma, grado_a));
      }
      a = a->get_next();
      b = b->get_next();
    }
    else if (grado_a < grado_b) { 
      nuevo = new SllPolyNode(pair_double_t(val_a, grado_a));
      a = a->get_next();
    }
    else {
      nuevo = new SllPolyNode(pair_double_t(val_b, grado_b));
      b = b->get_next();
    }

    if (nuevo != NULL) {
      if (sllpolsum.empty()) {
        sllpolsum.push_front(nuevo);
        tail = sllpolsum.get_head();
      } else {
        tail->set_next(nuevo);
        tail = nuevo;
      }
    }
  }

  while (a != NULL) {
    SllPolyNode* nuevo = new SllPolyNode(pair_double_t(
      a->get_data().get_val(),
      a->get_data().get_inx()
    ));

    if (sllpolsum.empty()) {
      sllpolsum.push_front(nuevo);
      tail = sllpolsum.get_head();
    } else {
      tail->set_next(nuevo);
      tail = nuevo;
    }

    a = a->get_next();
  }

  while (b != NULL) {
    SllPolyNode* nuevo = new SllPolyNode(pair_double_t(
      b->get_data().get_val(),
      b->get_data().get_inx()
    ));

    if (sllpolsum.empty()) {
      sllpolsum.push_front(nuevo);
      tail = sllpolsum.get_head();
    } else {
      tail->set_next(nuevo);
      tail = nuevo;
    }

    b = b->get_next();
  }
}

// MODIFICACION
double SllPolynomial::suma_coef() const {
  double suma{0.0};
  SllPolyNode* aux = get_head();
  while (aux != NULL) {
    suma += aux->get_data().get_val();
    aux = aux->get_next();
  }
  return suma;
}

#endif  // SLLPOLYNOMIAL_H_