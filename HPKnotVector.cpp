/**
 *  @file         HPKnotVector.cpp
 *  @brief        Class for managing a knot vector in the HPIGX project.
 *
 *  This class provides functionality for creating and manipulating
 *  a knot vector used in various numerical methods, particularly in
 *  spline and finite element representations. It includes methods
 *  for insertion, elevation, and retrieval of knot vector properties.
 *
 *  @author       Jin Lingzhi
 *  @email        jinlz0428@outlook.com
 *  @version      1.0.0
 *  @date         2024-09-01
 *  @copyright    2024 - 2027 Jin Lingzhi. All Rights Reserved.
 */

 // Include necessary header files
#include "HPKnotVector.h"

namespace hpigx {

	HPKnotVector::HPKnotVector(const VectorXd& knots, const index_t& degree)
		: m_value(knots.data(), knots.data() + knots.size()),
		m_degree(degree)
	{}

	HPKnotVector::HPKnotVector(const std::vector<value_t>& knots, const index_t& degree)
		: m_value(knots), m_degree(degree)
	{}

	const index_t HPKnotVector::sfind(const value_t& u) const {

		if (u < m_value[m_degree]) {
			return m_degree;
		}
		else if (u >= m_value[m_value.size() - m_degree]) {
			return m_value.size() - m_degree - 2;
		}

		auto it = std::upper_bound(m_value.begin() + m_degree, m_value.begin() + m_value.size() - m_degree, u);
		return std::distance(m_value.begin(), it) - 1;
	}

	const index_t HPKnotVector::ifind(const value_t& u) const {

		std::vector<value_t> uniqueKnots;
		std::unique_copy(m_value.begin(), m_value.end(), std::back_inserter(uniqueKnots));

		if (u < m_value.front())
			return 0;
		else if (u >= uniqueKnots.back())
			return uniqueKnots.size() - 2;

		auto it = std::upper_bound(uniqueKnots.begin(), uniqueKnots.end(), u);
		return std::distance(uniqueKnots.begin(), it) - 1;
	}

	const index_t HPKnotVector::efind(const value_t& u) const {

		std::vector<value_t> uniqueKnots;
		std::unique_copy(m_value.begin() + m_degree, m_value.end() - m_degree,
			std::back_inserter(uniqueKnots));

		if (u < m_value.front())
			return 0;
		else if (u >= uniqueKnots.back())
			return uniqueKnots.size() - 2;

		auto it = std::upper_bound(uniqueKnots.begin(), uniqueKnots.end(), u);
		return std::distance(uniqueKnots.begin(), it) - 1;
	}

	const VectorXd HPKnotVector::greville() const {

		VectorXd anchor = VectorXd::Zero(m_value.size() - m_degree - 1);

		for (size_t i = 0; i < anchor.rows(); ++i) {
			for (size_t j = 0; j < m_degree; ++j) {
				anchor(i) += m_value[i + j + 1];
			}
		}
		anchor /= m_degree;
		return anchor;
	}

	const value_t HPKnotVector::greville(const index_t& index) const {
		value_t anchor = 0;

		for (size_t j = 0; j < m_degree; ++j) {
			anchor += m_value[index + j + 1];
		}
		anchor /= m_degree;
		return anchor;
	}

	const VectorXd HPKnotVector::unique() const {
		std::vector<value_t> uniqueKnots;
		std::unique_copy(m_value.begin() + m_degree, m_value.end() - m_degree,
			std::back_inserter(uniqueKnots));

		return VectorXd::Map(uniqueKnots.data(), uniqueKnots.size());
	}

	const VectorXn HPKnotVector::multiplicities() const {

		std::vector<size_t> mults;
		mults.reserve(m_value.size());

		// Calculate the multiplicities of the knots
		for (size_t j = 0; j < m_value.size(); ++j) {

			size_t count = 1;
			while (j + 1 < m_value.size() && m_value[j] == m_value[j + 1]) {
				++count; // Increment multiplicity
				++j;     // Skip duplicate knots
			}
			mults.push_back(count);
		}
		return VectorXn::Map(mults.data(), mults.size());
	}

	const size_t HPKnotVector::multiplicities(const index_t& index) const {

		size_t mults = 0, t = 0;

		// Calculate the multiplicities of the knots
		for (size_t j = 0; j < m_value.size(); ++j) {

			size_t count = 1;
			// Check subsequent knots for duplicates and calculate multiplicity
			while (j + 1 < m_value.size() && m_value[j] == m_value[j + 1]) {
				++count; // Increment multiplicity
				++j;     // Skip duplicate knots
			}

			// If the current index matches the target index,
			// return the multiplicity
			if (index == t) {
				return count;
			}
			else {
				++t; // Move to the next index
			}
		}
		return mults;
	}

	const size_t HPKnotVector::multiplicity(const value_t& u) const {

		index_t span = this->sfind(u);
		size_t count = 0;

		// Check if the knot value at the span index is equal to u
		if (m_value[span] == u) {
			// Check for duplicate knots going backwards,
			// incrementing the multiplicity count
			count = 1;
			for (size_t j = span; j > 0 && m_value[j] == m_value[j - 1]; --j) {
				++count; // Increment multiplicity
			}
		}
		return count;
	}


	void HPKnotVector::insert(const value_t& u, const size_t& t) {
		m_value.reserve(m_value.size() + t);
		m_value.insert(m_value.begin() + this->sfind(u), u);
	}


	void HPKnotVector::insert(MatrixXd& coefs, const value_t& u, const size_t& t) {

		const size_t nrows = coefs.rows(); // Number of rows in the coefficient matrix
		index_t span = this->sfind(u);     // Find the span for the knot u

		// Shift coefficients to make space for new knot
		MatrixXd coefs_temp(nrows + t, coefs.cols());
		coefs_temp.topRows(span - m_degree + 1) = coefs.topRows(span - m_degree + 1);

		// Update coefficients for the new knot
		for (size_t j = 0; j < t; ++j) {
			coefs_temp.middleRows(span + 1, nrows + j - span) = coefs.middleRows(span, nrows + j - span);

			for (size_t j = span - m_degree + 1; j <= span; ++j) {
				value_t alpha = m_value[j + m_degree] - m_value[j];

				if (alpha != 0) {
					alpha = (u - m_value[j]) / alpha; // Calculate alpha
					coefs_temp.row(j) = alpha * coefs.row(j) + (1 - alpha) * coefs.row(j - 1);
				}
				else {
					coefs_temp.row(j) = coefs.row(j);
				}
			}

			// Insert the new knot into the knot vector
			m_value.insert(m_value.begin() + span + 1, u);
			coefs = coefs_temp; // Update the coefficients
			++span;             // Move to the next span
		}
	}

	void HPKnotVector::refine(const size_t& u)
	{
		if (u <= 0) {
			throw std::invalid_argument("The number of points to insert must be greater than zero.");
		}

		for (index_t count = 0; count < u; ++count) {
			size_t maxGapIndex = 0;
			value_t maxGap = 0.0;

			for (size_t i = 0; i < m_value.size() - 1; ++i) {
				value_t gap = m_value[i + 1] - m_value[i];
				if (gap > maxGap) {
					maxGap = gap;
					maxGapIndex = i;
				}
			}
			value_t newPoint = (m_value[maxGapIndex] + m_value[maxGapIndex + 1]) / 2.0;
			m_value.insert(m_value.begin() + maxGapIndex+1, newPoint);
		}
	}


	void HPKnotVector::refine(const VectorXd& u) {
		// Check if the input vector u is empty
		if (u.size() == 0) return;
		m_value.reserve(m_value.size() + u.size());

		// Insert elements from u into m_value
		m_value.insert(m_value.end(), u.begin(), u.end());
		std::sort(m_value.begin(), m_value.end());
	}

	void HPKnotVector::refine(MatrixXd& coefs, const VectorXd& u) {
		// Check if the input vector u is empty
		if (u.size() == 0) return;

		// Sort the input vector u
		VectorXd u_sorted; HPUtils::sort(u, u_sorted);

		const size_t nrows = coefs.rows(); // Number of rows in the coefficient matrix
		const size_t n = u_sorted.size();  // Number of cols in the coefficient matrix
		size_t first = this->sfind(u_sorted(0));
		size_t last = this->sfind(u_sorted(n - 1)) + 1;

		std::vector<value_t> value_temp(m_value);
		m_value.clear();
		m_value.resize(nrows + n + m_degree + 1);

		MatrixXd coefs_temp(coefs);
		coefs = MatrixXd::Zero(nrows + n, coefs_temp.cols());

		for (size_t i = 0; i < first - m_degree + 1; ++i) coefs.row(i) = coefs_temp.row(i);
		for (size_t i = last - 1; i < nrows; ++i) coefs.row(i + n) = coefs_temp.row(i);
		for (size_t i = 0; i < first + 1; ++i) m_value[i] = value_temp[i];
		

		for (size_t i = 0; i < value_temp.size(); ++i) {
			std::cout << value_temp[i];
			if (i < value_temp.size() - 1) {
				std::cout << ", ";
			}
		}

		for (size_t i = last + m_degree; i < nrows + m_degree + 1; ++i) m_value[i + n] = value_temp[i];

		size_t iold = last + m_degree - 1;     // Index for old coefficients
		size_t inew = last + m_degree + n - 1; // Index for new coefficients

		// Iterate backwards over the new knots
		for (index_t i = n - 1; i >= 0; --i) {

			// Move old knots into the new coefficient matrix
			while (u_sorted(i) <= value_temp[iold] && iold > first) {
				coefs.row(inew - m_degree - 1) = coefs_temp.row(iold - m_degree - 1);
				m_value[inew] = value_temp[iold];
				--iold;
				--inew;
			}
			coefs.row(inew - m_degree - 1) = coefs.row(inew - m_degree);

			// Update coefficients based on the new knot
			for (size_t j = 1; j <= m_degree; ++j) {
				size_t  index = inew - m_degree + j; // Calculate the index for the coefficients
				value_t alpha = m_value[inew + j] - u_sorted(i);

				if (alpha != 0) {
					alpha = alpha / (m_value[inew + j] - value_temp[iold - m_degree + j]); // Calculate alpha
					coefs.row(index - 1) = alpha * coefs.row(index - 1) + (1.0 - alpha) * coefs.row(index);
				}
				else {
					coefs.row(index - 1) = coefs.row(index);
				}
			}

			m_value[inew] = u_sorted(i); // Insert the new knot
			--inew;
		}
	}

	void HPKnotVector::refine(MatrixXd& coefs, const size_t& n)
	{
		std::vector<value_t> old_value = m_value;
		VectorXd diff = this->diff(n);
		m_value = old_value;
		this->refine(coefs, diff);
	}

	VectorXd HPKnotVector::diff(const size_t& n)
	{

		std::vector<value_t> old_value = m_value;
		this->refine(n);
		std::vector<double> differences;
		std::set_difference(
			m_value.begin(), m_value.end(),
			old_value.begin(), old_value.end(),
			std::back_inserter(differences));
		VectorXd diff = Eigen::Map<Eigen::VectorXd>(differences.data(), differences.size());
		return diff;
	}

	void HPKnotVector::remove(const value_t& u, const size_t& t) {

	}

	void HPKnotVector::remove(MatrixXd& coefs, const value_t& u,
		const size_t& t)
	{

	}

	void HPKnotVector::elevate(const size_t& t) {
		// Get unique knots from the knot vector
		VectorXd uniqueKnots = unique();
		m_value.reserve(m_value.size() + uniqueKnots.size() * t);

		// Insert each unique knot t times
		for (size_t j = 0; j < t; ++j) {
			m_value.insert(m_value.end(), uniqueKnots.begin(), uniqueKnots.end());
		}
		std::sort(m_value.begin(), m_value.end());
		m_degree += t; // Increase the degree by t
	}

	void HPKnotVector::elevate(MatrixXd& coefs, const size_t& t) {
		// Return if no elevation is needed
		if (t == 0) return;

		const size_t nrows = coefs.rows();
		const size_t ncols = coefs.cols();
		size_t p = m_degree;

		// Initialize control point matrices
		std::vector<MatrixXd> P(p + 1);
		for (size_t j = 0; j < p + 1; ++j) {
			P[j].setZero(nrows - j, ncols);
		}
		P[0] = coefs; // Set initial coefficients

		// Degree elevation using Cox-de Boor formula
		for (size_t i = 1; i <= p; ++i) {
			for (size_t j = 0; j < nrows - i; ++j) {
				if (m_value[j + p + 1] > m_value[i + j]) {
					P[i].row(j).noalias() =
						(P[i - 1].row(j + 1) - P[i - 1].row(j)) / (m_value[j + p + 1] - m_value[i + j]);
				}
			}
		}

		// Get knot multiplicities
		VectorXn mults = this->multiplicities();
		this->elevate(t); // Elevate knots

		const size_t nrows_new = m_value.size() - (p + t + 1);
		size_t p_new = p + t;

		// Initialize new control point matrices
		std::vector<MatrixXd> Q(p_new + 1);
		for (size_t i = 0; i < p_new + 1; ++i) {
			Q[i].setZero(nrows_new - i, ncols);
		}

		// Calculate scaling factors
		VectorXd factor = VectorXd::Ones(p + 1);
		for (size_t i = 0; i <= p; ++i) {
			for (size_t j = 1; j <= i; ++j) {
				factor[i] *= static_cast<value_t>(p + 1 - j) / static_cast<value_t>(p_new + 1 - j);
			}
		}

		index_t betak = 0; // sum of interior multiplicities
		for (size_t i = 0; i < mults.size() - 1; ++i) {
			// Set known coefficients
			for (size_t j = p + 1 - mults[i]; j <= p; ++j) {
				Q[j].row(betak + i * t) = factor[j] * P[j].row(betak);
			}

			for (index_t j = p_new - 1; j >= 0; --j) {
				// Fill triangular table
				for (size_t k = 1; k <= p_new - j; ++k) {
					// Update index k for the considered knot m_value
					size_t ik = k + betak + i * t;
					if (m_value[ik + p_new] > m_value[ik + j])
						Q[j].row(ik).noalias() =
						Q[j].row(ik - 1) + Q[j + 1].row(ik - 1) * (m_value[ik - 1 + p_new + 1] - m_value[ik - 1 + j + 1]);
				}
			}
			betak += mults(i + 1);
		}

		//Insert new coefficients into m_coefs
		coefs = MatrixXd::Map(Q[0].data(), Q[0].rows(), Q[0].cols()); // Map the data as modifiable
	}

	const size_t HPKnotVector::size() const {
		return m_value.size();
	}

	const VectorXd HPKnotVector::source() const {
		return VectorXd::Map(m_value.data(), m_value.size());
	}

	const size_t HPKnotVector::degree() const {
		return m_degree;
	}

	const size_t HPKnotVector::number() const {
		return m_value.size() - m_degree - 1;
	}

	const value_t HPKnotVector::minimum() const {
		return m_value[m_degree];
	}

	const value_t HPKnotVector::maximum() const {
		return m_value[m_value.size() - m_degree - 1];
	}

} // namespace hpigx