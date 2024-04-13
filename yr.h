/* 
 Yuksel Rating
 Author: temen
 Date(YYYY-MM-DD): 2024-04-13

 Description:
 This "library" helps perform skill-based matchmaking between two parties
 called "players". This is based off of Cem Yuksel's (cem@cemyuksel.com)
 "Skill-Based Matchmaking for Competitive Two-Player Games", ACM 2024.
 I am in no way associated with Cem Yuksel nor their works.

 See the LICENSE at the end of the file. 

 Before #including, 
     #define YUKSEL_RATING_IMPLEMENTATION 
 in the file / translation unit that you want to have the implementation. 
*/

#ifndef YUKSEL_RATING_H
#define YUKSEL_RATING_H

#define YUKSEL_RATING_DEFAULT_INITIAL_RATING 1500.0
#define YUKSEL_RATING_DEFAULT_ALPHA 2.0
#define YUKSEL_RATING_DEFAULT_MAX_RATING_CHANGE 350.0
#define YUKSEL_RATING_DEFAULT_RATING_CHANGE_SCALE 1.0
#define YUKSEL_RATING_DEFAULT_NUMBER_OF_GAMES_TRACKED 5
#define YUKSEL_RATING_DEFAULT_EXPECTED_WIN_RATE 0.5

#ifdef YUKSEL_RATING_STATIC
#define YUKSEL_RATING_DEF static
#else
#define YUKSEL_RATING_DEF extern
#endif

#ifdef __cplusplus
extern "C" {
#endif

/*
    _____   __________    __  ______  ______
   /  _/ | / / ____/ /   / / / / __ \/ ____/
   / //  |/ / /   / /   / / / / / / / __/   
 _/ // /|  / /___/ /___/ /_/ / /_/ / /___   
/___/_/ |_/\____/_____/\____/_____/_____/   
*/

typedef struct{
	double r; /* Rating */
	double W; /* Running weighted sum of variance weights */
	double R; /* Running weighted rating mean */
	double V; /* Running weighted rating variance */
	double D; /* Running weighted sum of derivatives of win probability */
	double w; /* Number of wins that are considered for tracking win rate */
	double m; /* Number of outcomes that are considered for tracking win rate */
}yuksel_player;

/* 
 Rationale for floating-point wins and outcomes: When computing the rating
 window, this assumes that it is a singular player that the window
 surrounds. In order to account for the fact that a "player" may represent
 a party of of players, you would set the wins and the outcomes to the an
 average of each of them respectively.
 Example: 5 players that have never queued before in party. 5 games tracked
     Ratings: 1250, 1375, 1500, 1625, 1750
     Wins:     2, 2, 3, 4, 5
     Outcomes: 3, 4, 4, 5, 5
     `player`: r = 1500, W = 0, R = 0, V = 0, D = 0, w = 3.2, m = 4.2 
*/

/*
 Initializes a player that has no game outcomes
 `a` is a single player that is being initialized
 `rating` is the initial rating value for the player
*/
YUKSEL_RATING_DEF void yuksel_rating_player_init(yuksel_player *a, double rating);

/*
 Initializes a player from a party of players
 `out` is a player that is being initialized
 `in` is the array of players within the party
 `player_count` is the number of players within the party
 `n` is the number of outcomes that are considered for tracking win rate
*/
YUKSEL_RATING_DEF void yuksel_rating_player_init_from_party(yuksel_player *out, yuksel_player *in, int player_count, int n);

/*
 Updates the rating of a player according to the score and their opponent
 `a` is the player that is being updated
     Assumed to be well-formed
 `b` is the opposing player `a`
     Assumed to be well-formed
 `s` is the score of player `a`
     Assumed to be is in [0, 1]
 `alpha` controls rating prediction responsivity to actual rating changes
     A good default: `alpha` = 2.0
 `delta_r_max` is the maximum change in rating in one update
     A good default: `delta_r_max` = 350.0
 `delta_r_scale` is the scale of the rating change
     A good default: `delta_r_scale` = 1.0
*/
YUKSEL_RATING_DEF void yuksel_rating_update(yuksel_player *a, const yuksel_player *b, double s, double alpha, double delta_r_max, double delta_r_scale);

/*
 Updates the rating of two players according to the score and each other
 `a` is a player opposing player `b` that is being updated
     Assumed to be well-formed
 `b` is a player opposing player `a` that is being updated
     Assumed to be well-formed
 `s` is the score of player `a`
     Assumed to be is in [0, 1]
 `alpha` controls rating prediction responsivity to actual rating changes
     A good default: `alpha` = 2.0
 `delta_r_max` is the maximum change in rating in one update
     A good default: `delta_r_max` = 350.0
 `delta_r_scale` is the scale of the rating change
     A good default: `delta_r_scale` = 1.0
*/
YUKSEL_RATING_DEF void yuksel_rating_update_zero_sum(yuksel_player *a, yuksel_player *b, double s, double alpha, double delta_r_max, double delta_r_scale);

/*
 Computes a rating window around a player to maintain win rate
 `r_star0` is the lower bound rating output of the rating window
 `r_star1` is the upper bound rating output of the rating window
 `a` is the player the rating window is being computed for
     Assumes `a->w` <= `n`, `a->m` <= `n`, and `a->w` <= `a->m`
 `n` is the number of outcomes that are considered for tracking win rate
 `lambda` is the desired win rate
     A good default: `lambda` = 0.5
*/
YUKSEL_RATING_DEF void yuksel_rating_rating_window(double *r_star0, double *r_star1, const yuksel_player *a, int n, double lambda);

#ifdef __cplusplus
}
#endif

#endif /* YUKSEL_RATING_H */

/*
    ____                __                          __        __  _           
   /  _/___ ___  ____  / /__  ____ ___  ___  ____  / /_____ _/ /_(_)___  ____ 
   / // __ `__ \/ __ \/ / _ \/ __ `__ \/ _ \/ __ \/ __/ __ `/ __/ / __ \/ __ \
 _/ // / / / / / /_/ / /  __/ / / / / /  __/ / / / /_/ /_/ / /_/ / /_/ / / / /
/___/_/ /_/ /_/ .___/_/\___/_/ /_/ /_/\___/_/ /_/\__/\__,_/\__/_/\____/_/ /_/ 
             /_/                                                              
*/

#ifdef YUKSEL_RATING_IMPLEMENTATION

#include <math.h>
#include <float.h>

void yuksel_rating_player_init(yuksel_player *a, double rating){
	a->r = rating;
	a->W = 0.0;
	a->R = 0.0;
	a->V = 0.0;
	a->D = 0.0;
	a->w = 0.0;
	a->m = 0.0;
}

void yuksel_rating_player_init_from_party(yuksel_player *out, yuksel_player *in, int player_count, int n){
	out->r = 0.0;
	out->W = 0.0;
	out->R = 0.0;
	out->V = 0.0;
	out->D = 0.0;
	out->w = 0.0;
	out->m = 0.0;

	if(player_count < 1 || n < 0)
		return;

	double N = (double)n;
	for(int i = 0; i < player_count; i++){
		out->r += in[i].r;
		out->W += in[i].W;
		out->R += in[i].R;
		out->V += in[i].V;
		out->D += in[i].D;
		double m = in[i].m > N ? N : in[i].m;
		double w = in[i].w > m ? m : in[i].w;
		out->w += w;
		out->m += m;
	}

	double count = (double)player_count;
	out->r /= count;
	out->W /= count;
	out->R /= count;
	out->V /= count;
	out->D /= count;
	out->w /= count;
	out->m /= count;
}

/* f_E(a, b) = 1 / (1 + 10^(-(r_a - r_b) / 400)) */
double yuksel_rating_elo_win_probability(double r_a, double r_b){
	return 1.0 / (1.0 + pow(10.0, -(r_a - r_b) / 400.0));
}

#if 0
/* d/dr_a f_E(a, b) = q f_E(a, b) (1 - f_E(a, b)) */
double yuksel_rating_elo_win_probability_derivative(double r_a, double r_b){
	const double q = 0.005756462732485115; /* ln(10) / 400 */
	double tmp = elo_win_probability(r_a, r_b);
	return q * tmp * (1.0 - tmp);
}
#endif

/* Similar to above but the parameter is the estimated win probability */
double yuksel_rating_elo_win_probability_derivative2(double p){
	const double q = 0.005756462732485115; /* ln(10) / 400 */
	return q * p * (1.0 - p);
}

#if 0
/* 1 / sqrt(1 + (3 q^2 phi^2) / pi^2) */
double yuksel_rating_glicko_rating_deviation_contribution_scale(double phi){
	const double q = 0.005756462732485115; /* ln(10) / 400 */
	return 1.0 / sqrt(1.0 + ((3.0 * q * q) * phi * phi) / (M_PI * M_PI));
}
#endif

/* Similar to above but the parameter is squared */
double yuksel_rating_glicko_rating_deviation_contribution_scale2(double phi_phi){
	const double q = 0.005756462732485115; /* ln(10) / 400 */
	return 1.0 / sqrt(1.0 + ((3.0 * q * q) * phi_phi) / (M_PI * M_PI));
}

#define YUKSEL_RATING_F yuksel_rating_elo_win_probability
#define YUKSEL_RATING_F_PRIME yuksel_rating_elo_win_probability_derivative2
#define YUKSEL_RATING_G yuksel_rating_glicko_rating_deviation_contribution_scale2

void yuksel_rating_update(yuksel_player *a, const yuksel_player *b, double s, double alpha, double delta_r_max, double delta_r_scale){
	/* Compute the arguments for the contribution scale */
	double phi_phi_b = b->W >= DBL_EPSILON ? b->V / b->W : 0.0;
	double phi_phi_a = a->W >= DBL_EPSILON ? a->V / a->W : 0.0;
	double alpha_alpha_phi_phi_a = alpha * alpha * phi_phi_a;

	/* The parameter is squared within these functions, pass them in squared */
	double g_of_phi_b = YUKSEL_RATING_G(phi_phi_b);
	double g_of_alpha_phi_a = YUKSEL_RATING_G(alpha_alpha_phi_phi_a);

	/* Evaluate f(a, b) and d/dr_a f(a, b) */
	double p_a = YUKSEL_RATING_F(a->r, b->r);
	double p_a_prime = YUKSEL_RATING_F_PRIME(p_a);

	/* Rating update */
	double F_a = g_of_phi_b * (s - p_a); /* Equ 17, not stored */
	double D_a = g_of_phi_b * p_a_prime + g_of_alpha_phi_a * a->D; /* Equ 15 */
	double delta_r = F_a / D_a; /* Not stored */
	if(delta_r_max >= DBL_EPSILON){
		if(delta_r > delta_r_max){
			delta_r = delta_r_max;
			D_a = F_a / delta_r;
		}
		else if(delta_r < -delta_r_max){
			delta_r = -delta_r_max;
			D_a = F_a / delta_r;
		}
	}
	a->D = D_a;

	/* Update stats */
	a->r = a->r + delta_r_scale * delta_r; /* Equ 6 */
	double omega = g_of_alpha_phi_a;
	a->W = omega * a->W + 1.0; /* Equ 9*/
	double delta_R = a->r - a->R;
	a->R = a->R + delta_R / a->W; /* Equ 10 */
	a->V = omega * a->V + delta_R * (a->r - a->R); /* Equ 11 */
}

void yuksel_rating_update_zero_sum(yuksel_player *a, yuksel_player *b, double s, double alpha, double delta_r_max, double delta_r_scale){
	/* Compute the arguments for the contribution scale */
	double phi_phi_b = b->W >= DBL_EPSILON ? b->V / b->W : 0.0;
	double phi_phi_a = a->W >= DBL_EPSILON ? a->V / a->W : 0.0;
	double alpha_alpha_phi_phi_a = alpha * alpha * phi_phi_a;
	double alpha_alpha_phi_phi_b = alpha * alpha * phi_phi_b;

	/* The parameter is squared within these functions, pass them in squared */
	double g_of_phi_a = YUKSEL_RATING_G(phi_phi_a);
	double g_of_phi_b = YUKSEL_RATING_G(phi_phi_b);
	double g_of_alpha_phi_a = YUKSEL_RATING_G(alpha_alpha_phi_phi_a);
	double g_of_alpha_phi_b = YUKSEL_RATING_G(alpha_alpha_phi_phi_b);
		
	/* Evaluate f(a, b) and d/dr_a f(a, b) */
	double p_a = YUKSEL_RATING_F(a->r, b->r);
	double p_a_prime = YUKSEL_RATING_F_PRIME(p_a);

	/* Rating update */
	double F_a = g_of_phi_b * (s - p_a); /* Equ 17, not stored */
	double F_b = g_of_phi_b * (p_a - s); /* Equ 17, not stored */
	double D_a = g_of_phi_b * p_a_prime + g_of_alpha_phi_a * a->D; /* Equ 15 */
	double D_b = g_of_phi_a * p_a_prime + g_of_alpha_phi_b * a->D; /* Equ 15 */
	double delta_r = (D_a * F_a - D_b * F_b) / (D_a * D_a + D_b * D_b);
	if(delta_r_max >= DBL_EPSILON){
		if(delta_r > delta_r_max){
			delta_r = delta_r_max;
			D_a = F_a / delta_r;
			D_b = -F_b / delta_r;
		}
		else if(delta_r < -delta_r_max){
			delta_r = -delta_r_max;
			D_a = F_a / delta_r;
			D_b = -F_b / delta_r;
		}
	}
	a->D = D_a;
	b->D = D_b;

	/* Update stats for player `a` */
	a->r = a->r + delta_r_scale * delta_r; /* Equ 6 */
	double omega = g_of_alpha_phi_a;
	a->W = omega * a->W + 1.0; /* Equ 9*/
	double delta_R = a->r - a->R;
	a->R = a->R + delta_R / a->W; /* Equ 10 */
	a->V = omega * a->V + delta_R * (a->r - a->R); /* Equ 11 */

	/* Update stats for player `b` */
	b->r = b->r + delta_r_scale * -delta_r; /* Equ 6 */
	omega = g_of_alpha_phi_a;
	b->W = omega * b->W + 1.0; /* Equ 9*/
	delta_R = b->r - b->R;
	b->R = b->R + delta_R / b->W; /* Equ 10 */
	b->V = omega * b->V + delta_R * (b->r - b->R); /* Equ 11 */
}

#undef YUKSEL_RATING_G
#undef YUKSEL_RATING_F_PRIME
#undef YUKSEL_RATING_F

/* r_b = r_a + 400 log_10(1 / p - 1) */
double yuksel_rating_elo_win_probability_to_rating(double p, double r_a){
	return r_a + 400.0 * log10(1.0 / p - 1.0);
}

void yuksel_rating_rating_window(double *r_star0, double *r_star1, const yuksel_player *a, int n, double lambda){
	if(n < 0){
		*r_star0 = a->r;
		*r_star1 = a->r;
		return;
	}
	double np1 = (double)(n + 1);
	double p = (lambda * (np1 + a->m) - a->w) / np1;
	double delta_p = 1.0 / (np1 + a->m);
	*r_star0 = yuksel_rating_elo_win_probability_to_rating(p - delta_p, a->r);
	*r_star1 = yuksel_rating_elo_win_probability_to_rating(p + delta_p, a->r);
}

#endif /* YUKSEL_RATING_IMPLEMENTATION */

/* 
 LICENSE

 This is free and unencumbered software released into the public domain.

 Anyone is free to copy, modify, publish, use, compile, sell, or
 distribute this software, either in source code form or as a compiled
 binary, for any purpose, commercial or non-commercial, and by any
 means.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR
 OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
 ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 OTHER DEALINGS IN THE SOFTWARE.
*/

