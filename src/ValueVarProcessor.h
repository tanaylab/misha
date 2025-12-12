#ifndef VALUEVARPROCESSOR_H_
#define VALUEVARPROCESSOR_H_

#include "TrackExpressionVars.h"
#include "GInterval.h"
#include "rdbutils.h"
#include <vector>

class ValueVarProcessor {
public:
	ValueVarProcessor(rdb::IntervUtils &iu) : m_iu(iu) {}

	void process_value_vars(
		TrackExpressionVars::Value_vars &value_vars,
		const GInterval &interval,
		unsigned idx
	);

private:
	[[maybe_unused]] rdb::IntervUtils &m_iu;

	void process_single_value_var(
		TrackExpressionVars::Value_var &var,
		const GInterval &interval,
		unsigned idx
	);
};

#endif /* VALUEVARPROCESSOR_H_ */

