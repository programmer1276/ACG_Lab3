# import numpy as np
# import matplotlib.pyplot as plt
# import tkinter as tk
# from tkinter import messagebox, filedialog
#
# def compute_illuminance(
#     W_mm, H_mm,
#     Wres, Hres,
#     xL, yL, zL,
#     I0,
#     cx, cy, R
# ):
#     # координатная сетка, центр области в (0,0)
#     xs = np.linspace(-W_mm / 2.0, W_mm / 2.0, Wres)
#     ys = np.linspace(-H_mm / 2.0, H_mm / 2.0, Hres)
#     X, Y = np.meshgrid(xs, ys)
#
#     # расстояние до источника
#     r2 = (X - xL) ** 2 + (Y - yL) ** 2 + zL ** 2
#
#     # ламбертовский источник:
#     # I(θ) = I0 cosθ
#     # E = I(θ) cosθ / r^2 = I0 cos^2θ / r^2
#     # cosθ = zL / r => E = I0 * zL^2 / r^4
#     E = I0 * (zL ** 2) / (r2 ** 2)
#
#     # нормировка в 0..255
#     E_max = E.max()
#     if E_max <= 0:
#         E_norm = np.zeros_like(E, dtype=np.uint8)
#     else:
#         E_norm = (E / E_max * 255.0).astype(np.uint8)
#
#     # маска круга
#     circle_mask = (X - cx) ** 2 + (Y - cy) ** 2 <= R ** 2
#     circle_values = E[circle_mask]
#     if circle_values.size > 0:
#         Emin_c = float(circle_values.min())
#         Emax_c = float(circle_values.max())
#         Eavg_c = float(circle_values.mean())
#     else:
#         Emin_c = Emax_c = Eavg_c = 0.0
#
#     # 5 контрольных точек
#     def point_E(x, y):
#         r2_pt = (x - xL) ** 2 + (y - yL) ** 2 + zL ** 2
#         return float(I0 * (zL ** 2) / (r2_pt ** 2))
#
#     control_points = {
#         "Центр круга": point_E(cx, cy),
#         "Пересечение с +X": point_E(cx + R, cy),
#         "Пересечение с -X": point_E(cx - R, cy),
#         "Пересечение с +Y": point_E(cx, cy + R),
#         "Пересечение с -Y": point_E(cx, cy - R),
#     }
#
#     return {
#         "E": E,
#         "E_norm": E_norm,
#         "xs": xs,
#         "ys": ys,
#         "circle_stats": {
#             "E_min": Emin_c,
#             "E_max": Emax_c,
#             "E_avg": Eavg_c,
#         },
#         "control_points": control_points,
#     }
#
# def on_compute():
#     try:
#         W_mm = float(entry_W.get())
#         H_mm = float(entry_H.get())
#         Wres = int(entry_Wres.get())
#         Hres = int(entry_Hres.get())
#         xL = float(entry_xL.get())
#         yL = float(entry_yL.get())
#         zL = float(entry_zL.get())
#         I0 = float(entry_I0.get())
#         cx = float(entry_cx.get())
#         cy = float(entry_cy.get())
#         R = float(entry_R.get())
#     except ValueError:
#         messagebox.showerror("Ошибка", "Проверьте ввод числовых параметров.")
#         return
#
#     if zL <= 0:
#         messagebox.showerror("Ошибка", "zL должен быть > 0 (источник над плоскостью).")
#         return
#
#     res = compute_illuminance(
#         W_mm, H_mm,
#         Wres, Hres,
#         xL, yL, zL,
#         I0,
#         cx, cy, R
#     )
#
#     global last_result
#     last_result = res
#
#     E = res["E"]
#     E_norm = res["E_norm"]
#
#     # изображение освещённости
#     plt.figure("Распределение освещенности (нормированное)")
#     plt.imshow(
#         E_norm,
#         cmap="gray",
#         origin="lower",
#         extent=[-W_mm / 2.0, W_mm / 2.0, -H_mm / 2.0, H_mm / 2.0],
#     )
#     plt.xlabel("X, мм")
#     plt.ylabel("Y, мм")
#     plt.title("Нормированное распределение освещенности")
#     plt.colorbar(label="Уровень (0-255)")
#     plt.tight_layout()
#
#     # горизонтальное сечение через центр области
#     center_row = E.shape[0] // 2
#     xs = res["xs"]
#     section = E[center_row, :]
#
#     plt.figure("Сечение через центр области")
#     plt.plot(xs, section)
#     plt.xlabel("X, мм")
#     plt.ylabel("E, усл. ед.")
#     plt.title("Освещенность вдоль горизонтального сечения через центр")
#     plt.grid(True)
#     plt.tight_layout()
#
#     # характеристики в messagebox + консоль
#     circle_stats = res["circle_stats"]
#     cp = res["control_points"]
#
#     lines = []
#     lines.append("Пять контрольных точек (E, усл. ед.):")
#     for name, val in cp.items():
#         lines.append(f"{name}: {val:.6e}")
#     lines.append("")
#     lines.append("Статистика по кругу:")
#     lines.append(f"E_min: {circle_stats['E_min']:.6e}")
#     lines.append(f"E_max: {circle_stats['E_max']:.6e}")
#     lines.append(f"E_avg: {circle_stats['E_avg']:.6e}")
#
#     txt = "\n".join(lines)
#     print(txt)
#     messagebox.showinfo("Результаты расчёта", txt)
#
#     plt.show()
#
# def on_save():
#     global last_result
#     if last_result is None:
#         messagebox.showwarning("Внимание", "Сначала выполните расчёт.")
#         return
#
#     E_norm = last_result["E_norm"]
#     file_path = filedialog.asksaveasfilename(
#         defaultextension=".png",
#         filetypes=[("PNG files", "*.png"), ("All files", "*.*")]
#     )
#     if not file_path:
#         return
#
#     plt.figure()
#     plt.imshow(E_norm, cmap="gray", origin="lower")
#     plt.axis("off")
#     plt.savefig(file_path, bbox_inches="tight", pad_inches=0)
#     plt.close()
#     messagebox.showinfo("Сохранение", f"Изображение сохранено:\n{file_path}")
#
# # ---- GUI ----
# root = tk.Tk()
# root.title("ЛР3: Расчет освещенности от точечного источника (Ламбертовская диаграмма)")
#
# # область
# tk.Label(root, text="Ширина области W, мм").grid(row=0, column=0, sticky="w")
# entry_W = tk.Entry(root)
# entry_W.grid(row=0, column=1)
# entry_W.insert(0, "2000")
#
# tk.Label(root, text="Высота области H, мм").grid(row=1, column=0, sticky="w")
# entry_H = tk.Entry(root)
# entry_H.grid(row=1, column=1)
# entry_H.insert(0, "2000")
#
# tk.Label(root, text="Разрешение Wres, пикс").grid(row=2, column=0, sticky="w")
# entry_Wres = tk.Entry(root)
# entry_Wres.grid(row=2, column=1)
# entry_Wres.insert(0, "400")
#
# tk.Label(root, text="Разрешение Hres, пикс").grid(row=3, column=0, sticky="w")
# entry_Hres = tk.Entry(root)
# entry_Hres.grid(row=3, column=1)
# entry_Hres.insert(0, "400")
#
# # источник
# tk.Label(root, text="xL, мм").grid(row=0, column=2, sticky="w")
# entry_xL = tk.Entry(root)
# entry_xL.grid(row=0, column=3)
# entry_xL.insert(0, "0")
#
# tk.Label(root, text="yL, мм").grid(row=1, column=2, sticky="w")
# entry_yL = tk.Entry(root)
# entry_yL.grid(row=1, column=3)
# entry_yL.insert(0, "0")
#
# tk.Label(root, text="zL, мм (>0)").grid(row=2, column=2, sticky="w")
# entry_zL = tk.Entry(root)
# entry_zL.grid(row=2, column=3)
# entry_zL.insert(0, "1000")
#
# tk.Label(root, text="I0, Вт/ср").grid(row=3, column=2, sticky="w")
# entry_I0 = tk.Entry(root)
# entry_I0.grid(row=3, column=3)
# entry_I0.insert(0, "1000")
#
# # круг
# tk.Label(root, text="x центра круга cx, мм").grid(row=4, column=0, sticky="w")
# entry_cx = tk.Entry(root)
# entry_cx.grid(row=4, column=1)
# entry_cx.insert(0, "0")
#
# tk.Label(root, text="y центра круга cy, мм").grid(row=5, column=0, sticky="w")
# entry_cy = tk.Entry(root)
# entry_cy.grid(row=5, column=1)
# entry_cy.insert(0, "0")
#
# tk.Label(root, text="Радиус круга R, мм").grid(row=6, column=0, sticky="w")
# entry_R = tk.Entry(root)
# entry_R.grid(row=6, column=1)
# entry_R.insert(0, "800")
#
# # кнопки
# btn_compute = tk.Button(root, text="Рассчитать и визуализировать", command=on_compute)
# btn_compute.grid(row=7, column=0, columnspan=2, pady=10, sticky="we")
#
# btn_save = tk.Button(root, text="Сохранить изображение", command=on_save)
# btn_save.grid(row=7, column=2, columnspan=2, pady=10, sticky="we")
#
# last_result = None
#
# if __name__ == "__main__":
#     root.mainloop()

from matplotlib.gridspec import GridSpec
import numpy as np
import tkinter as tk
from tkinter import messagebox, filedialog

import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from matplotlib.patches import Circle


def compute_illuminance(
    W_mm, H_mm,
    Wres, Hres,
    xL, yL, zL,
    I0,
    cx, cy, R
):
    # координатная сетка, центр области в (0,0)
    xs = np.linspace(-W_mm / 2.0, W_mm / 2.0, Wres)
    ys = np.linspace(-H_mm / 2.0, H_mm / 2.0, Hres)
    X, Y = np.meshgrid(xs, ys)

    # расстояние до источника
    r2 = (X - xL) ** 2 + (Y - yL) ** 2 + zL ** 2

    # ламбертовский источник:
    # I(θ) = I0 cosθ
    # E = I(θ) cosθ / r^2 = I0 cos^2θ / r^2
    # cosθ = zL / r => E = I0 * zL^2 / r^4
    E = I0 * (zL ** 2) / (r2 ** 2)

    # нормировка в 0..255 для визуализации
    E_max = E.max()
    if E_max <= 0:
        E_norm = np.zeros_like(E, dtype=np.uint8)
    else:
        E_norm = (E / E_max * 255.0).astype(np.uint8)

    # маска круга
    circle_mask = (X - cx) ** 2 + (Y - cy) ** 2 <= R ** 2
    circle_values = E[circle_mask]
    if circle_values.size > 0:
        Emin_c = float(circle_values.min())
        Emax_c = float(circle_values.max())
        Eavg_c = float(circle_values.mean())
    else:
        Emin_c = Emax_c = Eavg_c = 0.0

    # 5 контрольных точек
    def point_E(x, y):
        r2_pt = (x - xL) ** 2 + (y - yL) ** 2 + zL ** 2
        return float(I0 * (zL ** 2) / (r2_pt ** 2))

    control_points = {
        "Центр круга": point_E(cx, cy),
        "Пересечение с +X": point_E(cx + R, cy),
        "Пересечение с -X": point_E(cx - R, cy),
        "Пересечение с +Y": point_E(cx, cy + R),
        "Пересечение с -Y": point_E(cx, cy - R),
    }

    return {
        "E": E,
        "E_norm": E_norm,
        "xs": xs,
        "ys": ys,
        "circle_stats": {
            "E_min": Emin_c,
            "E_max": Emax_c,
            "E_avg": Eavg_c,
        },
        "control_points": control_points,
    }


# ----------------- GUI + встроенный matplotlib -----------------

root = tk.Tk()
root.title("ЛР3: Расчет освещенности от точечного источника (Ламбертовская диаграмма)")

# Основной фрейм: слева управление, справа графика
frame_controls = tk.Frame(root, padx=10, pady=10)
frame_controls.grid(row=0, column=0, sticky="n")

frame_plot = tk.Frame(root, padx=10, pady=10)
frame_plot.grid(row=0, column=1, sticky="nsew")

root.geometry("1400x900")
frame_plot.pack_propagate(False)


root.grid_columnconfigure(1, weight=1)
root.grid_rowconfigure(0, weight=1)

# ----- Поля ввода -----

# область
tk.Label(frame_controls, text="Ширина области W, мм").grid(row=0, column=0, sticky="w")
entry_W = tk.Entry(frame_controls, width=10)
entry_W.grid(row=0, column=1)
entry_W.insert(0, "2000")

tk.Label(frame_controls, text="Высота области H, мм").grid(row=1, column=0, sticky="w")
entry_H = tk.Entry(frame_controls, width=10)
entry_H.grid(row=1, column=1)
entry_H.insert(0, "2000")

tk.Label(frame_controls, text="Разрешение Wres, пикс").grid(row=2, column=0, sticky="w")
entry_Wres = tk.Entry(frame_controls, width=10)
entry_Wres.grid(row=2, column=1)
entry_Wres.insert(0, "400")

tk.Label(frame_controls, text="Разрешение Hres, пикс").grid(row=3, column=0, sticky="w")
entry_Hres = tk.Entry(frame_controls, width=10)
entry_Hres.grid(row=3, column=1)
entry_Hres.insert(0, "400")

# источник
tk.Label(frame_controls, text="xL, мм").grid(row=4, column=0, sticky="w")
entry_xL = tk.Entry(frame_controls, width=10)
entry_xL.grid(row=4, column=1)
entry_xL.insert(0, "0")

tk.Label(frame_controls, text="yL, мм").grid(row=5, column=0, sticky="w")
entry_yL = tk.Entry(frame_controls, width=10)
entry_yL.grid(row=5, column=1)
entry_yL.insert(0, "0")

tk.Label(frame_controls, text="zL, мм (>0)").grid(row=6, column=0, sticky="w")
entry_zL = tk.Entry(frame_controls, width=10)
entry_zL.grid(row=6, column=1)
entry_zL.insert(0, "1000")

tk.Label(frame_controls, text="I0, Вт/ср").grid(row=7, column=0, sticky="w")
entry_I0 = tk.Entry(frame_controls, width=10)
entry_I0.grid(row=7, column=1)
entry_I0.insert(0, "1000")

# круг
tk.Label(frame_controls, text="cx, мм (центр круга)").grid(row=8, column=0, sticky="w")
entry_cx = tk.Entry(frame_controls, width=10)
entry_cx.grid(row=8, column=1)
entry_cx.insert(0, "0")

tk.Label(frame_controls, text="cy, мм (центр круга)").grid(row=9, column=0, sticky="w")
entry_cy = tk.Entry(frame_controls, width=10)
entry_cy.grid(row=9, column=1)
entry_cy.insert(0, "0")

tk.Label(frame_controls, text="R, мм (радиус круга)").grid(row=10, column=0, sticky="w")
entry_R = tk.Entry(frame_controls, width=10)
entry_R.grid(row=10, column=1)
entry_R.insert(0, "800")

# Кнопки
def on_save():
    global last_result
    if last_result is None:
        messagebox.showwarning("Внимание", "Сначала выполните расчёт.")
        return

    E_norm = last_result["E_norm"]
    file_path = filedialog.asksaveasfilename(
        defaultextension=".png",
        filetypes=[("PNG files", "*.png"), ("All files", "*.*")]
    )
    if not file_path:
        return

    from matplotlib import pyplot as plt

    fig_save, ax = plt.subplots()
    ax.imshow(E_norm, cmap="gray", origin="lower")
    ax.axis("off")
    fig_save.savefig(file_path, bbox_inches="tight", pad_inches=0)
    plt.close(fig_save)
    messagebox.showinfo("Сохранение", f"Изображение сохранено:\n{file_path}")


btn_compute = tk.Button(frame_controls, text="Рассчитать и визуализировать")
btn_compute.grid(row=11, column=0, columnspan=2, pady=(10, 5), sticky="we")

btn_save = tk.Button(frame_controls, text="Сохранить изображение", command=on_save)
btn_save.grid(row=12, column=0, columnspan=2, pady=(0, 10), sticky="we")

# Поле для текстовых результатов
tk.Label(frame_controls, text="Результаты:").grid(row=13, column=0, columnspan=2, sticky="w")
results_text = tk.Text(frame_controls, width=40, height=12)
results_text.grid(row=14, column=0, columnspan=2, pady=(2, 0))


# ----- Фигура matplotlib внутри Tk -----

# fig = Figure(figsize=(7, 6), dpi=100)
#
# ax_img = fig.add_subplot(2, 1, 1)   # верх: карта освещённости
# ax_prof = fig.add_subplot(2, 1, 2)  # низ: сечение
fig = Figure(figsize=(12, 9), dpi=120)

# создаем сетку 2x1 с одинаковыми пропорциями строк
gs = GridSpec(2, 1, figure=fig, height_ratios=[1, 1])

ax_img = fig.add_subplot(gs[0])   # верхний — распределение
ax_prof = fig.add_subplot(gs[1])  # нижний — сечение

canvas = FigureCanvasTkAgg(fig, master=frame_plot)
canvas_widget = canvas.get_tk_widget()
canvas_widget.pack(fill="both", expand=True)

cbar = None  # для цветовой шкалы
last_result = None


# def update_plots(res, W_mm, H_mm, cx, cy, R):
#     global cbar
#
#     E = res["E"]
#     E_norm = res["E_norm"]
#     xs = res["xs"]
#
#     # ----- верхний график: карта освещенности -----
#     ax_img.clear()
#
#     extent = [-W_mm / 2.0, W_mm / 2.0, -H_mm / 2.0, H_mm / 2.0]
#
#     im = ax_img.imshow(
#         E_norm,
#         cmap="gray",
#         origin="lower",
#         extent=extent,
#         aspect="equal"
#     )
#
#     # круг области анализа
#     # circle = Circle((cx, cy), R, fill=False, linewidth=1.5)
#     # ax_img.add_patch(circle)
#
#     ax_img.set_title("Нормированное распределение освещенности")
#     ax_img.set_xlabel("X, мм")
#     ax_img.set_ylabel("Y, мм")
#
#     # цветовая шкала: пересоздаем при обновлении
#     if cbar is not None:
#         cbar.remove()
#     cbar = fig.colorbar(im, ax=ax_img, fraction=0.046, pad=0.04)
#     cbar.set_label("Уровень (0-255)")
#
#     # ----- нижний график: горизонтальное сечение -----
#     ax_prof.clear()
#     center_row = E.shape[0] // 2
#     section = E[center_row, :]
#
#     ax_prof.plot(xs, section)
#     ax_prof.set_xlabel("X, мм")
#     ax_prof.set_ylabel("E, усл. ед.")
#     ax_prof.set_title("Освещенность вдоль горизонтального сечения через центр")
#     ax_prof.grid(True)
#
#     canvas.draw()
def update_plots(res, W_mm, H_mm, cx, cy, R):
    global cbar

    E = res["E"]
    E_norm = res["E_norm"]
    xs = res["xs"]

    # ----- верхний график: карта освещенности -----
    ax_img.clear()

    extent = [-W_mm / 2.0, W_mm / 2.0, -H_mm / 2.0, H_mm / 2.0]

    im = ax_img.imshow(
        E_norm,
        cmap="gray",
        origin="lower",
        extent=extent,
        aspect="equal"
    )

    # если круг больше не нужен — ничего не добавляем
    # если хочешь вернуть круг:
    # from matplotlib.patches import Circle
    # circle = Circle((cx, cy), R, fill=False, linewidth=1.5, color='black')
    # ax_img.add_patch(circle)

    ax_img.set_title("Нормированное распределение освещенности")
    ax_img.set_xlabel("X, мм")
    ax_img.set_ylabel("Y, мм")

    # ---- цветовая шкала: создаем один раз, потом только обновляем ----
    if cbar is None:
        cbar = fig.colorbar(im, ax=ax_img, fraction=0.046, pad=0.04)
        cbar.set_label("Уровень (0-255)")
    else:
        # привязываем colorbar к новому изображению
        cbar.update_normal(im)

    # ----- нижний график: горизонтальное сечение -----
    ax_prof.clear()
    center_row = E.shape[0] // 2
    section = E[center_row, :]

    ax_prof.plot(xs, section)
    ax_prof.set_xlabel("X, мм")
    ax_prof.set_ylabel("E, усл. ед.")
    ax_prof.set_title("Освещенность вдоль горизонтального сечения через центр")
    ax_prof.grid(True)

    # перерисовка canvas
    canvas.draw()


def on_compute():
    try:
        W_mm = float(entry_W.get())
        H_mm = float(entry_H.get())
        Wres = int(entry_Wres.get())
        Hres = int(entry_Hres.get())
        xL = float(entry_xL.get())
        yL = float(entry_yL.get())
        zL = float(entry_zL.get())
        I0 = float(entry_I0.get())
        cx = float(entry_cx.get())
        cy = float(entry_cy.get())
        R = float(entry_R.get())
    except ValueError:
        messagebox.showerror("Ошибка", "Проверьте ввод числовых параметров.")
        return

    if zL <= 0:
        messagebox.showerror("Ошибка", "zL должен быть > 0 (источник над плоскостью).")
        return

    res = compute_illuminance(
        W_mm, H_mm,
        Wres, Hres,
        xL, yL, zL,
        I0,
        cx, cy, R
    )

    global last_result
    last_result = res

    # Обновляем графики
    update_plots(res, W_mm, H_mm, cx, cy, R)

    # Формируем текст результатов
    circle_stats = res["circle_stats"]
    cp = res["control_points"]

    lines = []
    lines.append("Пять контрольных точек (E, усл. ед.):")
    for name, val in cp.items():
        lines.append(f"{name}: {val:.6e}")
    lines.append("")
    lines.append("Статистика по кругу:")
    lines.append(f"E_min: {circle_stats['E_min']:.6e}")
    lines.append(f"E_max: {circle_stats['E_max']:.6e}")
    lines.append(f"E_avg: {circle_stats['E_avg']:.6e}")

    txt = "\n".join(lines)
    print(txt)

    results_text.delete("1.0", tk.END)
    results_text.insert(tk.END, txt)


# привязка обработчика к кнопке
btn_compute.config(command=on_compute)

if __name__ == "__main__":
    root.mainloop()
