#!/usr/bin/env julia

using Pkg

raw_selectors = isempty(ARGS) ?
    get(ENV, "ADAPTIVEOPTICS_TEST_SELECTORS", "") : join(ARGS, ',')
selectors = strip.(split(raw_selectors, ','))
any(isempty, selectors) && error("test selectors must be nonempty")
Pkg.test(test_args=selectors)
