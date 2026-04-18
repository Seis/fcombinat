#ifndef ABSYNTYPES
#define ABSYNTYPES
#include "node.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct Unit_s Unit;
typedef struct Id_s Id;
typedef struct Expr_s Expr;
typedef struct ExprList_s ExprList;
typedef struct Rule_s Rule;
typedef struct Rules_s Rules;
typedef struct Error_s Error;
typedef struct Spec_s Spec;

#include "parser.tab.h"

typedef enum { NONE, LESS, EQUAL, GREATER } Restriction;
typedef enum { LEXER, PARSER } ErrorType;
typedef enum { ISERROR, NOTERROR } SpecType;

struct Unit_s {
  enum yytokentype type;
  KeyT key;
  char *(*toString)(const struct Unit_s *self);
  char *(*toJson)(const struct Unit_s *self);
};

struct Id_s {
  char *name;
  KeyT key;
  char *(*toString)(const struct Id_s *self);
  char *(*toJson)(const struct Id_s *self);
};

struct Expr_s {
  void *component;
  enum yytokentype type;
  Restriction restriction;
  long long int limit;
  KeyT key;
  char *(*toString)(const struct Expr_s *self);
  char *(*toJson)(const struct Expr_s *self);
};

struct ExprList_s {
  Expr **components;
  int size;
  int space;
  KeyT key;
  char *(*toString)(const struct ExprList_s *self);
  char *(*toJson)(const struct ExprList_s *self);
};

struct Rule_s {
  Id *variable;
  Expr *expression;
  KeyT key;
  char *(*toString)(const struct Rule_s *self);
  char *(*toJson)(const struct Rule_s *self);
};

struct Rules_s {
  Rule **components;
  int size;
  int space;
  KeyT key;
  char *(*toString)(const struct Rules_s *self);
  char *(*toJson)(const struct Rules_s *self);
};

struct Error_s {
  int line;
  char *message;
  ErrorType type;
  KeyT key;
  char *(*toString)(const struct Error_s *self);
  char *(*toJson)(const struct Error_s *self);
};

struct Spec_s {
  SpecType type;
  void *component;
  KeyT key;
  char *(*toString)(const struct Spec_s *self);
  char *(*toJson)(const struct Spec_s *self);
};
#endif

#ifndef ABSYN_H
#define ABSYN_H

Unit *newUnit(enum yytokentype type);
Id *newId(char *name);
Expr *newExpr(void *component, enum yytokentype type, Restriction restriction,
              long long int limit);
ExprList *newExprList(Expr *expression);
Rule *newRule(Id *variable, Expr *expression);
Rules *newRules(Rule *statement);
Error *newError(int line, char *message, ErrorType type);
Spec *newSpec(void *component, SpecType type);

ExprList *addExprToList(Expr *expression, ExprList *list);
Rules *addRuleToList(Rule *statement, Rules *list);

void freeNodeRecursive(void *node, NodeType type);
void freeNode(void *node, NodeType type);

#endif