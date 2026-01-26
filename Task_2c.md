matrices_to_v first 'folds out' the four $2$ x $2$ matrices into four $8$-component vectors by calling on matrix_to_vector from task 2a for each matrix. The function transforms the matrices into vectors as stated in task 2a in the project description:

$$\begin{align*}

\begin{pmatrix}
    a_r+ia_i & b_r+ib_i \\
    c_r+ic_i & d_r+id_i
\end{pmatrix} \rightarrow \begin{pmatrix} a_r & b_r & c_r & d_r & a_i & b_i &c_i & d_i \end{pmatrix}
\\
\\
\begin{pmatrix}
    e_r+ie_i & f_r+if_i \\
    g_r+ig_i & h_r+ih_i
\end{pmatrix} \rightarrow \begin{pmatrix} e_r & f_r & g_r & h_r & e_i & f_i &g_i & h_i \end{pmatrix}
\\
\\
\begin{pmatrix}
    j_r+j_i & k_r+ik_i \\
    l_r+il_i & m_r+im_i
\end{pmatrix} \rightarrow \begin{pmatrix} j_r & k_r & l_r & m_r & j_i & k_i &l_i & m_i \end{pmatrix}
\\
\\
\begin{pmatrix}
    n_r+in_i & o_r+io_i \\
    p_r+ip_i & q_r+iq_i
\end{pmatrix} \rightarrow \begin{pmatrix} n_r & o_r & p_r & q_r & n_i & o_i &p_i & q_i \end{pmatrix}

\end{align*}$$

These four vectors are then transformed into a $32$-component vector $\vec{v}$ so that the real parts of the matrices are placed first, and the imaginary parts last:

$$\begin{align*}

\vec{v}=\begin{pmatrix} a_r \\ b_r \\ c_r \\ d_r \\ e_r \\ f_r \\ g_r \\ h_r \\ j_r \\ k_r \\ l_r \\ m_r \\ n_r \\ o_r \\ p_r \\ q_r \\ a_i \\ b_i \\ c_i \\ d_i \\ e_i \\ f_i \\ g_i \\ h_i \\ j_i \\ k_i \\ l_i \\ m_i \\ n_i \\ o_i \\ p_i \\ q_i \end{pmatrix}^T
\end{align*}$$