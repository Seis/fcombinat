/* 
	Parser for Combinatorial Spec described in http://www.maplesoft.com/support/help/maple/view.aspx?path=combstruct/specification
*/

%{
#include <stdio.h>
#include <string.h>
#include "grammar/absyn.h"

Spec* root; /* root of abstract syntax tree */
NodeST* ST; /* Symbol table with allocated nodes (for cleanup on parse error) */
int hasLexerError = 0;
int yyerror(char *msg);
extern int yylex();
extern FILE* yyin;
extern int lineNumber;
%}

%union 
{
  int symbol;
  Unit* unit;
  Id* id;
  int number;
  Expr* exp;
  ExprList* explist;
  Rule* stmt;
  Rules* stmtlist;
  Spec* grammar;
}

%token <unit> EPSILON 
%token <unit> ATOM
%token <unit> Z
%token <symbol> UNION
%token <symbol> PROD
%token <symbol> SET
%token <symbol> POWERSET
%token <symbol> SEQUENCE
%token <symbol> CYCLE
%token <symbol> SUBST
%token <symbol> CARD
%token <symbol> LPAR
%token <symbol> RPAR
%token <symbol> COMMA
%token <symbol> LEQ
%token <symbol> GEQ
%token <symbol> EQ
%token <symbol> LT
%token <symbol> GT
%token <symbol> LBRACK
%token <symbol> RBRACK
%token <id> ID
%token <number> NUMBER

%type <exp> expression 
%type <explist> expression_list
%type <stmt> statement
%type <stmtlist> statement_list
%type <grammar> grammar
%type <id> identifier

%start grammar

%%

grammar:	   		          statement_list { $$ = newSpec($1, NOTERROR) ; if (!hasLexerError) root = $$; }
;

statement_list:		                  statement { $$ = newRules($1); }
		 			| statement_list COMMA statement { $$ = addRuleToList($3, $1); }
;

statement:			          identifier EQ expression { $$ = newRule($1, $3); }
					| Z EQ expression { $$ = newRule(newId("Z"), $3); }
;

expression_list: 	                  expression { $$ = newExprList($1); }
		 			| expression_list COMMA expression { $$ = addExprToList($3, $1); }
;

identifier:               ID { $$ = $1; }
                        | ID LBRACK NUMBER RBRACK { 
                            char buffer[256];
                            snprintf(buffer, sizeof(buffer), "%s[%d]", $1->name, $3);
                            $$ = newId(buffer); 
                          }
;

expression:	   		          EPSILON { $$ = newExpr($1, EPSILON, NONE, 0); }
		 			| ATOM { $$ = newExpr($1, ATOM, NONE, 0); }
		 			| Z { $$ = newExpr($1, Z, NONE, 0); }
		 			| identifier { $$ = newExpr($1, ID, NONE, 0); }
		 			| UNION LPAR expression_list RPAR { $$ = newExpr($3, UNION, NONE, 0); }
                                        | PROD LPAR expression_list RPAR { $$ = newExpr($3, PROD, NONE, 0); }
		 			| SUBST LPAR expression_list RPAR { $$ = newExpr($3, SUBST, NONE, 0); }
		 			| SET LPAR expression RPAR { $$ = newExpr($3, SET, NONE, 0); }
                                        | SET LPAR expression COMMA CARD LEQ NUMBER RPAR { $$ = newExpr($3, SET, LESS, $7); }
		          	        | SET LPAR expression COMMA NUMBER GEQ CARD RPAR { $$ = newExpr($3, SET, LESS, $5); }
		 			| SET LPAR expression COMMA CARD EQ NUMBER RPAR { $$ = newExpr($3, SET, EQUAL, $7); }
		 			| SET LPAR expression COMMA NUMBER EQ CARD RPAR { $$ = newExpr($3, SET, EQUAL, $5); }
		 			| SET LPAR expression COMMA CARD GEQ NUMBER RPAR { $$ = newExpr($3, SET, GREATER, $7); }
		 			| SET LPAR expression COMMA NUMBER LEQ CARD RPAR { $$ = newExpr($3, SET, GREATER, $5); }
                    | SET LPAR expression COMMA CARD LT NUMBER RPAR { $$ = newExpr($3, SET, LESS, $7 - 1); }
                    | SET LPAR expression COMMA NUMBER GT CARD RPAR { $$ = newExpr($3, SET, LESS, $5 - 1); }
                    | SET LPAR expression COMMA CARD GT NUMBER RPAR { $$ = newExpr($3, SET, GREATER, $7 + 1); }
                    | SET LPAR expression COMMA NUMBER LT CARD RPAR { $$ = newExpr($3, SET, GREATER, $5 + 1); }
		 			| POWERSET LPAR expression RPAR { $$ = newExpr($3, POWERSET, NONE, 0); }
		 			| POWERSET LPAR expression COMMA CARD LEQ NUMBER RPAR { $$ = newExpr($3, POWERSET, LESS, $7); }
		 			| POWERSET LPAR expression COMMA NUMBER GEQ CARD RPAR { $$ = newExpr($3, POWERSET, LESS, $5); }
		 			| POWERSET LPAR expression COMMA CARD EQ NUMBER RPAR { $$ = newExpr($3, POWERSET, EQUAL, $7); }
		 			| POWERSET LPAR expression COMMA NUMBER EQ CARD RPAR { $$ = newExpr($3, POWERSET, EQUAL, $5); }
		 			| POWERSET LPAR expression COMMA CARD GEQ NUMBER RPAR { $$ = newExpr($3, POWERSET, GREATER, $7); }
		 			| POWERSET LPAR expression COMMA NUMBER LEQ CARD RPAR { $$ = newExpr($3, POWERSET, GREATER, $5); }
                    | POWERSET LPAR expression COMMA CARD LT NUMBER RPAR { $$ = newExpr($3, POWERSET, LESS, $7 - 1); }
                    | POWERSET LPAR expression COMMA NUMBER GT CARD RPAR { $$ = newExpr($3, POWERSET, LESS, $5 - 1); }
                    | POWERSET LPAR expression COMMA CARD GT NUMBER RPAR { $$ = newExpr($3, POWERSET, GREATER, $7 + 1); }
                    | POWERSET LPAR expression COMMA NUMBER LT CARD RPAR { $$ = newExpr($3, POWERSET, GREATER, $5 + 1); }
		 			| SEQUENCE LPAR expression RPAR { $$ = newExpr($3, SEQUENCE, NONE, 0); }
                 	                | SEQUENCE LPAR expression COMMA CARD LEQ NUMBER RPAR { $$ = newExpr($3, SEQUENCE, LESS, $7); }
		 			| SEQUENCE LPAR expression COMMA NUMBER GEQ CARD RPAR { $$ = newExpr($3, SEQUENCE, LESS, $5); }
                 	                | SEQUENCE LPAR expression COMMA CARD EQ NUMBER RPAR { $$ = newExpr($3, SEQUENCE, EQUAL, $7); }
		 			| SEQUENCE LPAR expression COMMA NUMBER EQ CARD RPAR { $$ = newExpr($3, SEQUENCE, EQUAL, $5); }
                 	                | SEQUENCE LPAR expression COMMA CARD GEQ NUMBER RPAR { $$ = newExpr($3, SEQUENCE, GREATER, $7); }
		 			| SEQUENCE LPAR expression COMMA NUMBER LEQ CARD RPAR { $$ = newExpr($3, SEQUENCE, GREATER, $5); }
                    | SEQUENCE LPAR expression COMMA CARD LT NUMBER RPAR { $$ = newExpr($3, SEQUENCE, LESS, $7 - 1); }
                    | SEQUENCE LPAR expression COMMA NUMBER GT CARD RPAR { $$ = newExpr($3, SEQUENCE, LESS, $5 - 1); }
                    | SEQUENCE LPAR expression COMMA CARD GT NUMBER RPAR { $$ = newExpr($3, SEQUENCE, GREATER, $7 + 1); }
                    | SEQUENCE LPAR expression COMMA NUMBER LT CARD RPAR { $$ = newExpr($3, SEQUENCE, GREATER, $5 + 1); }
                                        | CYCLE LPAR expression RPAR { $$ = newExpr($3, CYCLE, NONE, 0); }
                 	                | CYCLE LPAR expression COMMA CARD LEQ NUMBER RPAR { $$ = newExpr($3, CYCLE, LESS, $7); }
		 			| CYCLE LPAR expression COMMA NUMBER GEQ CARD RPAR { $$ = newExpr($3, CYCLE, LESS, $5); }
                 	                | CYCLE LPAR expression COMMA CARD EQ NUMBER RPAR { $$ = newExpr($3, CYCLE, EQUAL, $7); }
		 			| CYCLE LPAR expression COMMA NUMBER EQ CARD RPAR { $$ = newExpr($3, CYCLE, EQUAL, $5); }
                 	                | CYCLE LPAR expression COMMA CARD GEQ NUMBER RPAR { $$ = newExpr($3, CYCLE, GREATER, $7); }
		 			| CYCLE LPAR expression COMMA NUMBER LEQ CARD RPAR { $$ = newExpr($3, CYCLE, GREATER, $5); }
                    | CYCLE LPAR expression COMMA CARD LT NUMBER RPAR { $$ = newExpr($3, CYCLE, LESS, $7 - 1); }
                    | CYCLE LPAR expression COMMA NUMBER GT CARD RPAR { $$ = newExpr($3, CYCLE, LESS, $5 - 1); }
                    | CYCLE LPAR expression COMMA CARD GT NUMBER RPAR { $$ = newExpr($3, CYCLE, GREATER, $7 + 1); }
                    | CYCLE LPAR expression COMMA NUMBER LT CARD RPAR { $$ = newExpr($3, CYCLE, GREATER, $5 + 1); }
;


%%

int reportError(Error* error)
{
  if (error->type == LEXER) {
    hasLexerError = 1;
  }
  root = newSpec(error, ISERROR);
  char* str = error->toString(error);
  int result = fprintf(stderr, "%s\n", str);
  free(str);
  return result;
}

int yyerror(char *msg)
{
  return reportError(newError(lineNumber, msg, PARSER));
}

Spec* readSpec(char* filename)
{
  ST = newNodeST();
  yyin = fopen(filename, "r");
  yyparse();
  fclose(yyin);

  if (root->type == ISERROR) { // we need to clean up
  	// copy error info
  	Error* error = (Error*) root->component;
  	int line = error->line;
  	ErrorType type = error->type;
  	char* str = (char*) malloc(sizeof(char) * (strlen(error->message) + 1));
  	sprintf(str, "%s", error->message);

  	// free all nodes (including root) and their contents, and remove them from ST
  	cleanup(ST);

	// set root to error
  	root = newSpec(newError(line, str, type), ISERROR);
  }

  // free ST (but not abstract syntax tree nodes) since it is not needed anymore
  free(ST);

  return root;
}

// Flex definitions
typedef struct yy_buffer_state *YY_BUFFER_STATE;
extern YY_BUFFER_STATE yy_scan_string(const char *str);
extern void yy_delete_buffer(YY_BUFFER_STATE buffer);
extern void yy_switch_to_buffer(YY_BUFFER_STATE buffer);

Spec* readSpecStr(const char* input)
{
  ST = newNodeST();
  YY_BUFFER_STATE buffer = yy_scan_string(input);
  yy_switch_to_buffer(buffer);
  
  yyparse();
  
  yy_delete_buffer(buffer);

  if (root->type == ISERROR) { // we need to clean up
  	// copy error info
  	Error* error = (Error*) root->component;
  	int line = error->line;
  	ErrorType type = error->type;
  	char* str = (char*) malloc(sizeof(char) * (strlen(error->message) + 1));
  	sprintf(str, "%s", error->message);

  	// free all nodes (including root) and their contents, and remove them from ST
  	cleanup(ST);

	// set root to error
  	root = newSpec(newError(line, str, type), ISERROR);
  }

  // free ST (but not abstract syntax tree nodes) since it is not needed anymore
  free(ST);

  return root;
}

#ifndef _COMPILE_LIB
int main(int argc, char* argv[])
{
  readSpec(argv[1]);
  
  printf("%s\n", root->toJson(root));
}
#endif
