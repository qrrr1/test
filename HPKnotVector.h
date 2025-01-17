/**
 *  @file         HPKnotVector.h
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

#ifndef HPKNOTVECTOR_H
#define HPKNOTVECTOR_H

 // Include necessary header files
#include "HPPch.h"
#include "HPUtils.h"
#include "HPBaseIterator.h"

namespace hpigx {

    class HPKnotVector {

    public:
        // Type definitions for smart pointers
        using sptr = std::shared_ptr<HPKnotVector>;

        using uptr = std::unique_ptr<HPKnotVector>;

    public:
        /**
         * @brief Default copy constructor.
         *
         * @param other The other HPKnotVector to copy from.
         */
        HPKnotVector(const HPKnotVector&) = default;

        /**
         * @brief Default constructor.
         */
        HPKnotVector() = default;

        /**
         * @brief Constructs a knot vector with given knots and degree.
         *
         * @param knots The vector of knot values.
         * @param degree The degree of the knot vector.
         */
        HPKnotVector(const VectorXd& knots, const index_t& degree);

        HPKnotVector(const std::vector<value_t>& knots, const index_t& degree);

        /**
         * @brief Default destructor.
         */
        ~HPKnotVector() = default;

        /**
         * @brief Finds the span index for a given parameter value.
         *
         * @param u The parameter value to find the span for.
         * @return The span index for parameter `u`.
         */
        const index_t sfind(const value_t& u) const;

        /**
         * @brief Finds the interval index for a given parameter value.
         *
         * @param u The parameter value to find the knot for.
         * @return The interval index for parameter `u`.
         */
        const index_t ifind(const value_t& u) const;

        /**
         * @brief Finds the element index for a given parameter value.
         *
         * @param u The parameter value to find the end knot for.
         * @return The element index for parameter `u`.
         */
        const index_t efind(const value_t& u) const;

        /**
         * @brief Computes and returns the Greville abscissae.
         *
         * @return A matrix containing the Greville abscissae.
         */
        const VectorXd greville() const;

        /**
         * @brief Computes the Greville abscissa at a specific index.
         *
         * @param index The index of the Greville abscissa to compute.
         * @return The Greville abscissa at the specified index.
         */
        const value_t greville(const index_t& index) const;

        /**
         * @brief Returns a vector of unique knots.
         *
         * @return A vector of unique knot values.
         */
        const VectorXd unique() const;

        /**
         * @brief Returns the multiplicities of knots.
         *
         * @return A vector containing the multiplicities of each knot.
         */
        const VectorXn multiplicities() const;

        /**
         * @brief Returns the multiplicity of the knot at a specific index.
         *
         * @param index The index of the knot.
         * @return The multiplicity of the knot at `index`.
         */
        const size_t multiplicities(const index_t& index) const;

        /**
         * @brief Returns the multiplicity of a specific parameter value.
         *
         * @param u The parameter value to find the multiplicity for.
         * @return The multiplicity of the knot corresponding to parameter `u`.
         */
        const size_t multiplicity(const value_t& u) const;

        /**
         * @brief Returns the multiplicity of a specific parameter value.
         *
         * @param u The parameter value to find the multiplicity for.
         * @return The multiplicity of the knot corresponding to parameter `u`.
         */

    public:
        /**
         * @brief Inserts a new knot into the knot vector.
         *
         * @param u The parameter value of the knot to insert.
         * @param t The index at which to insert the knot.
         */
        void insert(const value_t& u, const size_t& t);

        /**
         * @brief Inserts a new knot and associated coefficients.
         *
         * @param coefs The matrix of coefficients to associate with the knot.
         * @param u The parameter value of the knot to insert.
         * @param t The index at which to insert the knot.
         */
        void insert(MatrixXd& coefs, const value_t& u, const size_t& t);

        /**
         * @brief Refines the knot vector by adding new knots at specified parameter values.
         *
         * @param u A vector of parameter values to refine the knot vector.
         */
        void refine(const VectorXd& u);

        void refine(const size_t& n);




        /**
         * @brief Refines the knot vector and updates associated coefficients.
         *
         * @param coefs The matrix of coefficients to update.
         * @param u A vector of parameter values to refine the knot vector.
         */
        void refine(MatrixXd& coefs, const VectorXd& u);

        void refine(MatrixXd& coefs, const size_t& n);

        VectorXd diff(const size_t& n);

        /**
         * @brief Removes a knot from the knot vector.
         *
         * @param u The parameter value of the knot to remove.
         * @param t The index at which to remove the knot.
         */
        void remove(const value_t& u, const size_t& t);

        /**
         * @brief Elevates the degree of the knot vector.
         *
         * @param t The number of degrees to elevate.
         */
        void remove(MatrixXd& coefs, const value_t& u, const size_t& t);

        /**
         * @brief Elevates the degree of the knot vector.
         *
         * @param t The number of degrees to elevate.
         */
        void elevate(const size_t& t);

        /**
         * @brief Elevates the degree of the knot vector and updates associated coefficients.
         *
         * @param coefs The matrix of coefficients to update.
         * @param t The number of degrees to elevate.
         */
        void elevate(MatrixXd& coefs, const size_t& t);

    public:
        /**
         * @brief Returns the number of knots in the knot vector.
         *
         * @return The size of the knot vector.
         */
        const size_t size() const;

        /**
         * @brief Returns the knot values as a vector.
         *
         * @return A vector containing the knot values.
         */
        const VectorXd source() const;

        /**
         * @brief Returns the degree of the knot vector.
         *
         * @return The degree of the knot vector.
         */
        const size_t degree() const;

        /**
         * @brief Returns the number of unique knots in the knot vector.
         *
         * @return The number of unique knots.
         */
        const size_t number() const;

        const value_t minimum() const;

        const value_t maximum() const;

    private:
        // Custom iterator class for the knot vector
        class Iterator : public HPBaseIterator<value_t> {

        private:
            typename std::vector<value_t>::iterator m_iter; // Using std::vector iterator

            index_t m_index;                                // Current index of the iterator
            size_t m_length;                                // Total length of the vector

        public:
            /**
             * @brief Constructs an iterator for the knot vector.
             *
             * @param iterator The underlying vector iterator.
             * @param index The current index of the iterator.
             * @param length The total length of the vector.
             */
            Iterator(typename std::vector<value_t>::iterator iterator, index_t index, size_t length)
                : m_iter(iterator), m_index(index), m_length(length) {}

            ~Iterator() = default;

            value_t operator*() override {
                return *m_iter; // Dereference the iterator
            }

            void operator++() override {
                if (isEnd()) {
                    return; // Already at the end, no need to increment
                }
                ++m_index;  // Increment index
                ++m_iter;   // Prefix increment
            }

            void operator++(int) override {
                if (isEnd()) {
                    return; // Already at the end, no need to increment
                }
                ++m_index;  // Increment index
                ++m_iter;   // Prefix increment
            }

            void operator--() override {
                if (isBegin()) {
                    return; // Already at the beginning, no need to decrement
                }
                --m_index;  // Decrement index
                --m_iter;   // Prefix decrement
            }

            bool operator!=(const HPBaseIterator<value_t>& other) const override {
                const Iterator* otherIt = static_cast<const Iterator*>(&other);
                if (!otherIt) {
                    throw std::invalid_argument("Comparison with incompatible iterator."); // Check for valid comparison
                }
                return m_iter != otherIt->m_iter; // Not equal comparison
            }

            bool operator==(const HPBaseIterator<value_t>& other) const override {
                const Iterator* otherIt = static_cast<const Iterator*>(&other);
                if (!otherIt) {
                    throw std::invalid_argument("Comparison with incompatible iterator."); // Check for valid comparison
                }
                return m_iter == otherIt->m_iter; // Equal comparison
            }

            value_t current() override {
                return *m_iter; // Return the current element
            }

            bool next() override {
                if (isEnd()) {
                    return false; // Cannot move to next if at the end
                }
                ++m_index;   // Increment index
                ++m_iter;    // Move to the next element
                return true; // Successfully moved to the next element
            }

            bool isBegin() const override {
                return m_index <= 0; // At the beginning if index is 0
            }

            bool isEnd() const override {
                return m_index >= static_cast<index_t>(m_length); // At the end if index is equal to or greater than length
            }
        };

    public:
        /**
         * @brief Returns the beginning iterator for the knot vector.
         *
         * @return An iterator pointing to the beginning of the knot vector.
         */
        typename HPKnotVector::Iterator begin() {
            return Iterator(m_value.begin(), 0, m_value.size());
        }

        /**
         * @brief Returns the end iterator for the knot vector.
         *
         * @return An iterator pointing to the end of the knot vector.
         */
        typename HPKnotVector::Iterator end() {
            return Iterator(m_value.end(), m_value.size(), m_value.size());
        }

        value_t& at(const index_t i) {
            return m_value[i];
        }

        const value_t& at(const index_t i) const {
            return m_value[i];
        }

        /**
         * @brief Overloads the assignment operator to assign knot values from a VectorXd.
         *
         * @param knots The vector of knot values to assign.
         * @return A reference to this HPKnotVector.
         */
        HPKnotVector& operator=(const VectorXd& knots) {
            m_value = std::vector<value_t>(knots.data(), knots.data() + knots.size());
            return *this;
        }

        /**
         * @brief Overloads the subscript operator to access knot values.
         *
         * @param i The index of the knot to access.
         * @return Reference to the knot value at index `i`.
         */
        value_t& operator()(index_t i) {
            return m_value[i];
        }

        /**
         * @brief Overloads the subscript operator to access knot values (const).
         *
         * @param i The index of the knot to access.
         * @return Reference to the knot value at index `i`.
         */
        const value_t& operator()(index_t i) const {
            return m_value[i];
        }

    protected:
        // Declaration of member variables
        std::vector<value_t> m_value; // Vector storing the knot values
        size_t m_degree;              // Degree of the knot vector
    };

} // namespace hpigx

#endif // HPKNOTVECTOR_H