#ifndef INTERVVARPROCESSOR_H_
#define INTERVVARPROCESSOR_H_

#include "TrackExpressionVars.h"
#include "GInterval.h"
#include "rdbutils.h"
#include <vector>

class IntervVarProcessor {
public:
	IntervVarProcessor(rdb::IntervUtils &iu) : m_iu(iu) {}

	void process_interv_vars(
		TrackExpressionVars::Interv_vars &interv_vars,
		const GInterval &interval,
		unsigned idx
	);

private:
	rdb::IntervUtils &m_iu;

	void process_distance(
		TrackExpressionVars::Interv_var &var,
		const GInterval &interval,
		unsigned idx
	);

	void process_distance_center(
		TrackExpressionVars::Interv_var &var,
		const GInterval &interval,
		unsigned idx
	);

	void process_distance_edge(
		TrackExpressionVars::Interv_var &var,
		const GInterval &interval,
		unsigned idx
	);

	void process_coverage(
		TrackExpressionVars::Interv_var &var,
		const GInterval &interval,
		unsigned idx
	);

	void process_neighbor_count(
		TrackExpressionVars::Interv_var &var,
		const GInterval &interval,
		unsigned idx
	);
};

#endif /* INTERVVARPROCESSOR_H_ */

