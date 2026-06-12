#pragma once

using namespace std;

void cal_def_grad_bar (vector<array<array<double, 3>, 3>>& delta_def_grad,
                       const vector<double>& det_delta_def_grad);