/* wr_metric.h */

/* WOHLGEMUTH-RONACHER SPIKE DISTANCE
 *
 * This code computes the spike distance metric described in 
 *     S. Wohlgemuth & B. Ronacher (2007), Auditory Discrimination of
 *     Amplitude Modulations Based On Metric Distances of Spike Trains,
 *     J. Neurophysiol. 97:3082-3092.*
 *
 * To compile this code type
 *
 *     mex -v wr_metric_mex.c wr_metric.c
 *
 * This implementation:
 * Author:    Robert Mill 
 * Version:   0.1 (12-12-12)
 * Revisions:
 *   <none>
 */

#ifndef WR_METRIC_H_INCLUDED
#define WR_METRIC_H_INCLUDED


/* MACRO DEFINITIONS */

#define MAX_ERR_LEN 128


/* TYPE DEFINITIONS */

typedef struct WR_SpikeTrains WR_SpikeTrains;


/* STRUCTURE DEFINITIONS */

/* Structure to hold pointers to spike train data */
struct WR_SpikeTrains
{
    /* Number of spike trains being stored */
    int n_trains;

    /* Spike counts for each spike train */
    int *np_counts;

    /* Spike times for each spike train */
    double **dpp_times_ms;
};


/* FUNCTION PROTOTYPES */

/**
 * Clear any error messages that have been set.
 */
void WR_ClearError( void );


/**
 * Set an error message. Calling this function will overwrite any error
 * message that was previously set.
 *
 * Parameters:
 *   cp_msg - pointer to null-terminated string
 */
void WR_SetError( const char *cp_msg );


/**
 * Get the error message. This function returns a pointer to the error
 * message that was most recently stored. If no error message has been
 * set (or it has been cleared), NULL is returned.
 *
 * Return:
 *   pointer to null-terminated string (the error message)
 *   or NULL if no error
 */
const char * WR_GetError( void );


/**
 * Check whether an error has occurred.
 *
 * Return:
 *   1 if error message is set
 *   0 otherwise
 */
int WR_IsError( void );


/**
 * Allocate memory for spike train (pointer) data. This function allocates
 * memory to store pointers to spike train data. Note that the function
 * does NOT allocate space to store the spike times. Rather it is assumed
 * that there are several buffers of doubles which contain spike times.
 * This function reserves enough space to store:
 *  - spike counts for each train (in np_counts)
 *  - pointers to the spike times themselves (in dpp_times_ms)
 *
 * If the routine fails, NULL is returned.
 *
 * Return:
 *    pointer to spike trains structure
 *    or NULL on failure
 */
WR_SpikeTrains * WR_AllocSpikeTrains( int n_trains );


/**
 * Free memory allocated to spike train (pointer) data. This function
 * should be used to free the memory allocated by calls to
 * WR_AllocSpikeTrains. It is safe to call this function with a null
 * argument.
 *
 * Parameters:
 *    sp_trains - pointer to spike trains structure
 */
void WR_FreeSpikeTrains( WR_SpikeTrains *sp_trains );


/**
 * Run a routine to compute the pairwise distances between spike trains.
 * Input should be supplied in the form of a spike trains structure
 * (see WR_AllocSpikeTrains). Output is written to a double-precision
 * matrix pointed to by dp_out. If there are N spike trains in sp_trains
 * then dp_out should point to a buffer with enough space to store N x N
 * spike train distances. It is the responsibility of the caller to
 * allocate this matrix.
 *
 * Parameters:
 *    sp_trains - pointer to a spike trains structure
 *    dp_out    - pointer to results buffer
 *    d_alpha   - "alpha" parameter
 */
void WR_ComputeDistance(
    const WR_SpikeTrains *sp_trains,
    double *dp_out,
    double d_alpha );

#endif

/* end of file */
