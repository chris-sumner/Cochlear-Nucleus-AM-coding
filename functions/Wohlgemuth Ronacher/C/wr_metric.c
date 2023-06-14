/* wr_metric.c */

/* WOHLGEMUTH-RONACHER DISTANCE METRIC
 * 
 * See wr_metric.h for details.
 */


/* INCLUDE HEADERS */

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "wr_metric.h"


/* GLOBAL VARIABLES */

static char cp_errormsg [MAX_ERR_LEN];


/* FUNCTION DEFINITIONS */

/* For procedure description see wr_metric.h */
void WR_ClearError( void )
{
    /* Clear error message */
    memset( cp_errormsg, (char) 0, MAX_ERR_LEN );
}


/* For procedure description see wr_metric.h */
void WR_SetError( const char *cp_msg )
{
    /* Set the error message */
    strncpy( cp_errormsg, cp_msg, MAX_ERR_LEN - 1 );
}


/* For procedure description see wr_metric.h */
const char * WR_GetError( void )
{
    /* Get the error message */
    return cp_errormsg;
}


/* For procedure description see wr_metric.h */
int WR_IsError( void )
{
    /* Test whether an error message has been set */
    return *cp_errormsg != 0;
}


/* For procedure description see wr_metric.h */
WR_SpikeTrains * WR_AllocSpikeTrains( int n_trains )
{
    /* Allocate space for the spike trains */
    WR_SpikeTrains *sp_trains = calloc( 1, sizeof( WR_SpikeTrains ) );

    /* Allocated structure OK? */
    if( sp_trains != NULL )
    {
        /* Set number of spike trains */
        sp_trains->n_trains = n_trains;

        /* Allocate memory for spike counts and times */
        sp_trains->np_counts    = calloc( n_trains, sizeof( int ) );
        sp_trains->dpp_times_ms = calloc( n_trains, sizeof( double * ) );

        /* Failed to allocate memory for all fields? */
        if( sp_trains->np_counts == NULL ||
            sp_trains->dpp_times_ms == NULL )
        {
            /* Free any memory that was allocated */
            free( sp_trains->np_counts );
            free( sp_trains->dpp_times_ms );
            free( sp_trains );
        }
    }

    /* Return a pointer to newly-allocated structure (or NULL) */
    return sp_trains;
}


/* For procedure description see wr_metric.h */
void WR_FreeSpikeTrains( WR_SpikeTrains *sp_trains )
{
    /* Spike train structure allocated? */
    if( sp_trains )
    {
        /* Free its fields */
        free( sp_trains->np_counts );
        free( sp_trains->dpp_times_ms );

        /* Free the memory allocated to the structure */
        free( sp_trains );
    }
}


/* For procedure description see wr_metric.h */
void WR_ComputeDistance(
    const WR_SpikeTrains *sp_trains,
    double *dp_out,
    double d_alpha )
{
    /* Local variables */
    int i, j;
    int n_trains;

    double d_now  = -DBL_MAX;
    double d_next = +DBL_MAX;

    /* Pointers to dynamically-allocated memory */
    int    *np_pos = NULL;
    double *dp_xs  = NULL;
    double *dp_ys  = NULL;

    /* At start (i.e. prior to) or end of spike train? */
    int n_infstart = 1;
    int n_infend   = 0;

    /* Precomputed values */
    double d_a, d_a2;
    double d_T, d_T2, d_2aT;
    double d_exp, d_exp2;

    /* Precompute alpha^2 */
    d_a  = d_alpha;
    d_a2 = d_a * d_a;

    /* Set number of spike trains locally */
    n_trains = sp_trains->n_trains;

    /* Dynamically allocate intermediate variables */
    np_pos = calloc( n_trains*n_trains, sizeof( int ) );
    dp_xs  = calloc( n_trains, sizeof( double ) );
    dp_ys  = calloc( n_trains, sizeof( double ) );

    /* Failed to allocate any memory? */
    if( np_pos == NULL || dp_xs == NULL || dp_ys == NULL )
    {
        /* Free all memory that was allocated */
        free( np_pos );
        free( dp_xs );
        free( dp_ys );

        /* Set error message and return */
        WR_SetError( "Cannot allocate memory for filter data.\n" );
        return;
    }

    /* For each unique spike time: */
    do
    {
        /* Assume next spike "infinitely" in the future */
        d_next   = DBL_MAX;
        n_infend = 1;

        /* For each spike train i: */
        for( i = 0; i < n_trains; i++ )
        {
            /* Are there any more spikes left in this train? */
            if( np_pos[i] < sp_trains->np_counts[i] )
            {
                /* Get the time of the next spike */
                double *dp_times_ms = sp_trains->dpp_times_ms[i];
                double  d_time_ms   = dp_times_ms[np_pos[i]];

                /* Is the next spike sooner than the soonest found? */
                if( d_time_ms > d_now && d_time_ms < d_next )
                {
                    /* Update the soonest future spike time */
                    d_next   = d_time_ms;
                    n_infend = 0;
                }
            }
        }

        /* Integrating over infinite time? */
        if( n_infstart || n_infend )
        {
            /* Set time of integration to zero and exponential to zero.
             * [This prevents NaN errors from 0 * Inf.]
             */
            d_T = d_T2 = d_2aT = d_exp = d_exp2 = 0;
        }
        else
        {
            /* Get elapsed time, precompute quantities used in integral */
            d_T    = d_next - d_now;
            d_T2   = d_T * d_T;
            d_2aT  = 2 * d_a * d_T;            
            d_exp  = exp( -d_a * d_T );
            d_exp2 = d_exp * d_exp;
        }

        /* For each spike train i: */
        for( i = 0; i < n_trains; i++ )
        {
            /* Block-local variables for integrals */
            double d_I1, d_I2, d_I3;

            /* Get state of filter for train i */
            double d_ki1 = dp_ys[i];
            double d_ki2 = dp_xs[i];

            /* Solve three integrals */
            d_I1 =  1/( 2 * d_a ) * ( 1 - d_exp2 );
            d_I2 =  1/( 4 * d_a2) * ( 1 - d_exp2 * ( 1 + d_2aT ) );
            d_I3 = -1/( 2 * d_a ) * d_exp2 * d_T*d_T + d_I2/d_a;

            /* Square, integrate, accumulate (see wr_metric.h) */
            dp_out[i + i*n_trains] +=
                d_ki1 * d_ki1 * d_I1 +
                d_ki1 * d_ki2 * d_I2 * 2 +
                d_ki2 * d_ki2 * d_I3; 

            /* For each other spike train j: */
            for( j = i + 1; j < n_trains; j++ )
            {
                /* Get state of filter for train j */
                double d_kj1 = dp_ys[j];
                double d_kj2 = dp_xs[j];

                /* Multiply, integrate, accumulate (see wr_metric.h) */
                dp_out[i + j*n_trains] +=
                    d_ki1 * d_kj1 * d_I1 +
                    d_ki2 * d_kj2 * d_I3 + 
                    ( d_ki1 * d_kj2 + d_kj1 * d_ki2 ) * d_I2;
            }
        }

        /* Move current time to spike of next spike */
        d_now = d_next;

        /* For each spike train i: */
        for( i = 0; i < n_trains; i++ )
        {
            /* Update filters */
            dp_ys[i] = d_exp * ( dp_ys[i] + dp_xs[i]*d_T );
            dp_xs[i] = d_exp * dp_xs[i];

            /* Are there any more spikes left in this train? */
            if( np_pos[i] < sp_trains->np_counts[i] )
            {
                /* Get the time of the next spike */
                double *dp_times_ms = sp_trains->dpp_times_ms[i];
                double  d_time_ms   = dp_times_ms[np_pos[i]];

                /* Is the next spike "now"? */
                if( d_time_ms == d_now )
                {
                    /* Increment filter and advance to next spike */
                    np_pos[i]++;
                    dp_xs[i]++;
                }
            }
        }

        /* Performed one iteration - previous spike not in infinite past */
        n_infstart = 0;
    }
    while( n_infend != 1 );

    /* For each spike train i: */
    for( i = 0; i < n_trains; i++ )
    {
        /* Get index in memory of diagonal component for i */
        int ix_i = i*n_trains + i;

        /* For each spike train j: */
        for( j = i + 1; j < n_trains; j++ )
        {
            /* Get index in memory of diagonal component for j */
            int ix_j = j + j*n_trains;

            /* Get index in memory of off-diagonal component for (i, j) */
            int ix_ij = i + j*n_trains;
            int ix_ji = j + i*n_trains;

            /* Compute square difference */
            dp_out[ix_ij] = dp_out[ix_i] + dp_out[ix_j] - 2*dp_out[ix_ij];
            dp_out[ix_ji] = dp_out[ix_ij];
        }

        /* Store zero along diagonal */
        dp_out[ix_i] = 0;
    }

    /* Free dynamically-allocated variables */
    free( np_pos );
    free( dp_xs );
    free( dp_ys );
}

/* end of file */
