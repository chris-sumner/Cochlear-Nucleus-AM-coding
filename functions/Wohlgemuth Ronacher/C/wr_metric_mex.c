/* wr_metric_mex.c */

/* MEX INTERFACE FOR THE WOHLGEMUTH-RONACHER SPIKE DISTANCE
 *
 * This code computes the spike distance metric described in 
 *     S. Wohlgemuth & B. Ronacher (2007), Auditory Discrimination of
 *     Amplitude Modulations Based On Metric Distances of Spike Trains,
 *     J. Neurophysiol. 97:3082-3092.*
 *
 * For implementation details see wr_metric.h.
 *
 * To compile this code type
 *
 *     mex -v wr_metric_mex.c wr_metric.c
 *
 * This implementation:
 * Author:    Robert Mill 
 * Version:   0.2 (17-12-12)
 * Revisions:
 *   17/12/12 RWM Made vector direction lenient. (Can be columns or rows.)
 */

/* INCLUDE HEADERS */

#include <mex.h>
#include "wr_metric.h"


/* FUNCTION DEFINITIONS */

/**
 * Exit the Mex routine with an error message. The error message is
 * retrieved using WR_GetError.
 */
void ExitWithError( void )
{
    /* Get error message and printed to Matlab console with exit */
    mexErrMsgTxt( WR_GetError( ) );
}


/**
 * Convert a Matlab cell array to a C pointer structure. This function
 * takes a pointer to a Matlab cell array containing spike train data.
 * Each cell should contain a column vector of double-precision spike
 * times. Any allocation, type or dimension errors cause the routine to
 * return NULL.
 *
 * Parameters:
 *    mp_spikecell - pointer to a spike train cell array
 * Return:
 *    a pointer to a C pointer structure (WR_SpikeTrains)
 *    or NULL on failure
 * Error:
 *    NULL is returned on error and the error message is set
 */
WR_SpikeTrains * GetSpikeTrains( const mxArray *mp_spikecell )
{
    int i;
    int n_trains = mxGetM( mp_spikecell );

    /* Allocate space for spike train structure */
    WR_SpikeTrains *sp_trains = WR_AllocSpikeTrains( n_trains );

    /* Not able to allocate spike train pointer structure? */
    if( sp_trains == NULL )
    {
        /* Generate an error message and return NULL */
        WR_SetError( "Could not allocate memory." );
        return NULL;
    }

    /* For each spike train i: */
    for( i = 0; i < n_trains; i++ )
    {
        /* Get the cell containing the spike train */
        mxArray *mp_train = mxGetCell( mp_spikecell, i );

        /* Is input some format other than double? */
        if( !mxIsClass( mp_train, "double" ) )
        {
            /* Free memory, set error message and return */
            WR_FreeSpikeTrains( sp_trains );
            WR_SetError( "Spike times should be double format." );
            return NULL;
        }

        /* Is the spike train not a vector? */
        if( mxGetM( mp_train ) > 1 && mxGetN( mp_train ) > 1 )
        {
            /* Free memory, set error message and return */
            WR_FreeSpikeTrains( sp_trains );
            WR_SetError( "Each spike train should be a vector." );
            return NULL;
        }

        /* Set pointer to data and number of elements */
        sp_trains->np_counts[i]    = mxGetNumberOfElements( mp_train );
        sp_trains->dpp_times_ms[i] = (double *) mxGetData( mp_train );
    }

    /* Return pointer to spike train structure */
    return sp_trains;
}


/**
 * Main entry point from Matlab into C. For interpretation of the arguments
 * see wr_dist.m and inline comments.
 */
void mexFunction(
    int n_lhs,
    mxArray *mp_lhs [],
    int n_rhs,
    const mxArray *mp_rhs [] )
{
    /* Local variables (input) */
    WR_SpikeTrains *sp_trains = NULL;
    double          d_alpha   = 0;

    /* Local variables (output) */
    mxArray *mp_out    = NULL;
    double  *dp_out    = NULL;

    /* Clear errors */
    WR_ClearError( );

    /* Failed to supply at least one input? */
    if( n_rhs < 2 )
    {
        /* Generate error message */
        WR_SetError( "This function requires two input arguments." );
        ExitWithError( );
    }

    /* Is input some format other than cell? */
    if( !mxIsClass( mp_rhs[0], "cell" ) )
    {
        /* Generate error message */
        WR_SetError( "First argument requires cell format." );
        ExitWithError( );
    }

    /* Is the input not a column? */
    if( mxGetM( mp_rhs[0] ) == 0 || mxGetN( mp_rhs[0] ) != 1 )
    {
        /* Generate error message */
        WR_SetError( "First argument should be a column." );
        ExitWithError( );
    }

    /* Is alpha not a double? */
    if( !mxIsClass( mp_rhs[1], "double" ) )
    {
        /* Generate error message */
        WR_SetError( "Second argument requires double format." );
        ExitWithError( );
    }

    /* Is alpha not a scalar? */
    if( mxGetNumberOfElements( mp_rhs[1] ) != 1 )
    {
        /* Generate error message */
        WR_SetError( "Second argument should be a scalar." );
        ExitWithError( );
    }

    /* Put spike trains into a standard C format */
    sp_trains = GetSpikeTrains( mp_rhs[0] );

    /* Errors when attempting this? */
    if( WR_IsError( ) )
    {
        /* Generate error message */
        ExitWithError( );
    }

    /* Get alpha scalar */
    d_alpha = mxGetScalar( mp_rhs[1] );

    /* Is alpha not positive? */
    if( d_alpha <= 0 )
    {
        /* Generate error message */
        WR_SetError( "Second argument should be (strictly) positive." );
        ExitWithError( );
    }

    /* Reserve space for output */
    mp_out = mxCreateDoubleMatrix(
        sp_trains->n_trains,
        sp_trains->n_trains, 
        mxREAL );

    /* Get pointer to output space */
    dp_out = (double *) mxGetData( mp_out );

    /* Compute distances between all pairs of spike trains */
    WR_ComputeDistance( sp_trains, dp_out, d_alpha );

    /* Free up the spike train pointer structures */
    WR_FreeSpikeTrains( sp_trains );

    /* Return pointer to output */
    mp_lhs[0] = mp_out;
}

/* end of file */
